#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#define MAXSTR 51200

typedef struct vcf_info vcf_info;
typedef struct contig_info contig_info;
typedef struct mut_info mut_info;
typedef struct gtf_info gtf_info;
typedef struct gtf_metadata gtf_metadata;

struct vcf_info {
        char chrom[20];
        int pos;
        char id[20];
        char ref[512];
        char alt[512];
        char qual[10];
        char filter[10];
        char info[512];
};
struct mut_info {
	char mutation[512];
	int pos;
	int offset;
	int flag;
};
struct contig_info {
	char id[50];
	int len;
	int mut_cnt;
	int *deleted_pos;
	int del_cnt;
	mut_info *mut;
};
struct gtf_metadata {
	char seqname[50];
	int cnt;
};
struct gtf_info {
	char seqname[50];
	char source[50];
	char feature[50];
	int st;
	int en;
	char score[10];
	char strand;
	char frame;
	char group[5120];
	char transcript_id[100];
	int dont_print;
};

int count_alt_vcf(FILE *, int *);
vcf_info *read_alt_vcf(FILE *, contig_info **, int, int);
char *move_ptr(char *,char);
void create_mut_table(vcf_info *, int, contig_info *, int);
void add_deleted_pos_memory(int **, int, int);
void add_deleted_transcript(char ***, int, char *); 
int check_gtf(FILE *, int *, gtf_metadata **, int *);
char **read_headers(FILE *, int);
void alter_description(char ***, int, char *);
void alter_group(char *, mut_info *, int);
void write_alt_gtf(FILE *, gtf_info *, int, contig_info);
void write_norm_gtf(FILE *, gtf_info *, int);
void parse_transcript_id(gtf_info *);
void reset_flags(mut_info *, int);
void read_and_write_gtf(FILE *, FILE *, contig_info *, int, gtf_metadata *, int);
void print_usage();
/* alt_liftover.c: performs 'liftover' for a gtf file which refer to a standard reference genome, 
so that it instead contains coordinates referencing an alternate genome (previously produced by 
GATK's FastaAlternateReferenceMaker). 
## NOTE: USE CURATED VCF FILE (produced by curate_vcf.c) so that you can be sure what 
FastaAlternateReferenceMaker has done! ## */
/************************************************************************************************/
void main(int argc,char **argv)
{
	FILE *f, *g, *h;
	contig_info *contigs = NULL;	
	vcf_info *vcf_lines = NULL;
	gtf_metadata *gtf = NULL;
	mut_info **mut = NULL;
	int meta_cnt = 0, header_cnt = 0, contig_cnt = 0, vcf_lines_cnt = 0, gtf_lines_cnt, i, j, c; 
	char stem[100], output[100], *pnt, **header = NULL, suffix[512], vcf_name[512], gtf_name[512];

	suffix[0] = vcf_name[0] = gtf_name[0] = '\0';
	while ((c = getopt (argc, argv, "s:g:v:h")) != -1) {
                switch (c) {
                        case 's':
                                strcpy(suffix, optarg);
                                break;
                        case 'g':
                                strcpy(gtf_name, optarg);
                                break;
                        case 'v':
                                strcpy(vcf_name, optarg);
                                break;
			case 'h':
				print_usage();
				exit(0);
                        default: /*'?'*/
				
				print_usage();
                                exit(0);
                }
        }
	if (suffix[0] == '\0') 
		strcpy(suffix, "_alt");
	if (vcf_name[0] == '\0' || gtf_name[0] == '\0') {
		print_usage();
		exit(0);
	}

	if ((f = fopen(vcf_name,"r")) == NULL) {
                printf("Can't open file %s\n",vcf_name);
                exit(0);
        }

	vcf_lines_cnt = count_alt_vcf(f, &contig_cnt);
        printf("Found %d contig lines and %d lines of mutation information\n", contig_cnt, vcf_lines_cnt);
	vcf_lines = read_alt_vcf(f, &contigs, vcf_lines_cnt, contig_cnt);
	for (i=0; i<contig_cnt; ++i)
		if ((contigs[i].mut=calloc(contigs[i].mut_cnt, sizeof(mut_info)))==NULL) {
         		printf("memory allocation error in main\n");
                	exit(0);
        	}

	create_mut_table(vcf_lines, vcf_lines_cnt, contigs, contig_cnt);
	fclose(f);
	if (vcf_lines != NULL)
		free(vcf_lines);

	for (i=0; i<contig_cnt; ++i) {
		printf("%d) %s %d %d: del_cnt %d\ndeletions: ", i, contigs[i].id, contigs[i].len, contigs[i].mut_cnt, contigs[i].del_cnt);
		for (j=0; j<contigs[i].del_cnt; ++j) 
			printf("%d,", contigs[i].deleted_pos[j]);
		printf("\n");
		for (j=0; j<contigs[i].mut_cnt; j++)
			printf("%d.%d) %s, %d, %d\n", i, j, contigs[i].mut[j].mutation, contigs[i].mut[j].pos, contigs[i].mut[j].offset);
	}
	if ((g = fopen(gtf_name,"r")) == NULL) {
                printf("Can't open file %s\n",gtf_name);
                exit(0);
        }
	gtf_lines_cnt = check_gtf(g, &header_cnt, &gtf, &meta_cnt);
        printf("Found %d header lines and %d lines of gtf information\n", header_cnt, gtf_lines_cnt);	
	header = read_headers(g, header_cnt);
	for (i = 0; i<meta_cnt; ++i) {
		printf("%d) %s %d\n", i, gtf[i].seqname, gtf[i].cnt);
	}	

        strcpy(stem, gtf_name);
        if ((pnt = strstr(stem, ".gtf")) != NULL)
                *pnt = '\0';
        else {
                printf("error: specify .gtf file - use .gtf extension\n");
		exit(0);
	}
	sprintf(output, "%s%s.gtf", stem, suffix);
        if ((h = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n",output);
                exit(0);
        }
	if (header_cnt > 0)
		alter_description(&header, header_cnt, vcf_name);
	for (i = 0; i<header_cnt; ++i) 
		if (header[i] != NULL)
			fprintf(h, "%s", header[i]);
	for (i=0; i<header_cnt; ++i) {
                if (header[i] != NULL)
                        free(header[i]);
        }
        if (header != NULL)
                free(header);
	
	read_and_write_gtf(g, h, contigs, contig_cnt, gtf, meta_cnt);

	if (gtf !=NULL)
		free(gtf);
	for (i = 0; i<contig_cnt; ++i) {
		if (contigs[i].deleted_pos != NULL)
			free(contigs[i].deleted_pos);
		if (contigs[i].mut != NULL)
			free(contigs[i].mut);
	}
	if (contigs != NULL)
		free(contigs);
	fclose(g);
	fclose(h);
	exit(0);
}
/************************************************************************************************/

/************************************************************************************************/
void print_usage() 
{
	printf("Usage:\nRequired\n\t-g\tSpecify gtf file for liftover\n");
	printf("\t-v\tSpecify curated vcf file used to produce the alternate genome\n");
	printf("Optional\n\t-s\tSpecify suffix for output gtf file (default is _alt)\n");
	printf("\t-h\tPrint help\n");
	return;
}
/************************************************************************************************/

/************************************************************************************************/
int count_alt_vcf(FILE *f, int *contig_cnt)
{
        int cnt = 0;
        char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#') {
			if (strstr(line, "contig=") != NULL)
                        	*contig_cnt+=1;
		}
                else 
                        cnt++;
        }
        rewind(f);

        return(cnt);
}
/************************************************************************************************/

/************************************************************************************************/
int check_gtf(FILE *f, int *header_cnt, gtf_metadata **gtf, int *meta_cnt)
{
        int cnt = 0, i;
        char line[MAXSTR], tmp[20]; 
	gtf_metadata *tmp_realloc = NULL;

	*meta_cnt = 0;
	tmp[0] = '\0';
	if ((*gtf = calloc(1, sizeof(gtf_metadata))) == NULL) {
		printf("Memory error in check_gtf\n");
		exit(0);
	}

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#')
                        *header_cnt+=1;
                else {
			tmp[0] = '\0';
			sscanf(line, "%[^\t]", tmp);
			if (cnt == 0) {
				strcpy((*gtf)[*meta_cnt].seqname, tmp);
				(*gtf)[*meta_cnt].cnt = 1;
				*meta_cnt+=1;
			} 	
			else {
				if (strcmp(tmp, (*gtf)[(*meta_cnt)-1].seqname) != 0) {
					for (i = 0; i<*meta_cnt; ++i) {
						if (strcmp(tmp, (*gtf)[i].seqname) == 0) {
							printf("Problem with gtf file: entries for chromosomes are not grouped together\n");
							exit(0);	
						}
					}
					*meta_cnt+=1;
					if ((tmp_realloc = realloc((*gtf), ((*meta_cnt) * sizeof(gtf_metadata)))) == NULL) {
                        			printf("memory allocation error in check_gtf() - realloc\n");
                        			exit(0);
                			}
                			else {
                        			*gtf = tmp_realloc;
					}	
                                	strcpy((*gtf)[(*meta_cnt)-1].seqname, tmp);
					(*gtf)[(*meta_cnt)-1].cnt = 0;
				}
				(*gtf)[(*meta_cnt)-1].cnt++;		
			}	
                        cnt++;
		}
        }
        rewind(f);
	printf("Check gtf file: %d contigs, all records grouped for each\n", *meta_cnt);
        return(cnt);
}
/************************************************************************************************/

/************************************************************************************************/
vcf_info *read_alt_vcf(FILE *f, contig_info **contigs, int lines_cnt, int contig_cnt)
{
        vcf_info *tmp = NULL;
        int i = 0, j = 0, k=0;
        char line[MAXSTR], *pnt = NULL, *ptr = NULL;

        if ((*contigs=calloc(contig_cnt, sizeof(contig_info)))==NULL) {
                printf("memory allocation error in read_alt_vcf\n");
                exit(0);
        }
	for (k = 0; k<contig_cnt; ++k)
		(*contigs)[k].mut_cnt = 0;
        
	if ((tmp=calloc(lines_cnt, sizeof(vcf_info)))==NULL) {
                printf("memory allocation error in read_alt_vcf\n");
                exit(0);
        }

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#') {
			if ((strstr(line, "contig=") != NULL) && i < contig_cnt) {
                        	ptr = strstr(line, "ID=");
				pnt = move_ptr(ptr, ',');
				strcpy((*contigs)[i].id, ptr + 3);
				ptr = pnt + 1;
				*pnt = ',';
				pnt = move_ptr(ptr, '>');
				(*contigs)[i].len = atoi(ptr + 7);
				*pnt = '>';
                        	i++;
                	}
		}
                else if (j < lines_cnt) {
                        sscanf(line, "%s %d %s %[^\t] %[^\t] %s %s %s", tmp[j].chrom, &tmp[j].pos, tmp[j].id, tmp[j].ref, tmp[j].alt, tmp[j].qual, tmp[j].filter, tmp[j].info);
			for (k = 0; k<contig_cnt; ++k) {
				if (strcmp(tmp[j].chrom, (*contigs)[k].id) == 0)
					(*contigs)[k].mut_cnt++;
			}
                        j++;
                }
        }
        return tmp;
}
/************************************************************************************************/

/************************************************************************************************/
char **read_headers(FILE *f, int cnt)
{
	int i = 0;
	char **tmp = NULL;
        char line[MAXSTR];

        if ((tmp=calloc(cnt, sizeof(char *)))==NULL) {
                printf("memory allocation error in read_header\n");
                exit(0);
        }
        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#' && i < cnt) {
			if ((tmp[i]=calloc((strlen(line) + 1), sizeof(char)))==NULL) {
                                printf("memory allocation error in read_header\n");
                                exit(0);
                        }
                        strcpy(tmp[i], line);
                        i++;
		}
	}
	rewind(f);
	return tmp;
}
/************************************************************************************************/

/************************************************************************************************/
char *move_ptr(char *str,char delim)
{
char *tmp = NULL;

if((tmp = strchr(str,delim)) != NULL) *tmp = '\0';

return(tmp);
}
/************************************************************************************************/

/************************************************************************************************/
/*
	'pos' has to be ref genome pos - same as v.pos
	'offset' is new, calculated based on mutation - v.alt and v.ref
	'mutation' is info based on alt genome
	save deleted position numbers (based on ref genome) in 'contigs'
*/
void create_mut_table(vcf_info *v, int lines_cnt, contig_info *contigs, int contig_cnt)
{
	int i, j, k = 0, l, offset = 0, tmp_offset = 0, *tmp_delpos;
	
	for (i=0; i<contig_cnt; ++i) {
		offset=0;
		tmp_offset = 0;
		for (j = 0; j<contigs[i].mut_cnt; ++j) {
			contigs[i].mut[j].pos = v[k].pos;
			sprintf(contigs[i].mut[j].mutation, "%s->%s@%d", v[k].ref, v[k].alt, (v[k].pos + offset));
			tmp_offset = strlen(v[k].alt) - strlen(v[k].ref);
			offset += tmp_offset;
			if (tmp_offset < 0) {
				tmp_offset = -tmp_offset;
				add_deleted_pos_memory(&(contigs[i].deleted_pos), contigs[i].del_cnt, tmp_offset);
				for (l=0; l<tmp_offset; l++) {
					contigs[i].deleted_pos[(contigs[i].del_cnt + l)] = v[k].pos + 1 + l;
				} 
				contigs[i].del_cnt += tmp_offset;
			}
			contigs[i].mut[j].offset = offset;
			k++;
		}
	}
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void add_deleted_pos_memory(int **del, int old, int new) 
{
	int *tmp = NULL;

        if (*del == NULL) {
                if ((*del = calloc(new, sizeof(int))) == NULL) {
                        printf("memory allocation error in add_deleted_pos_memory()\n");
                        exit(0);
                }
        }
        else {
                if ((tmp = realloc(*del, ((old + new) * sizeof(int)))) == NULL) {
                        printf("memory allocation error in add_deleted_pos_memory() - realloc\n");
                        exit(0);
                }
                else
                        *del = tmp;
        }
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void add_deleted_transcript(char ***del, int old, char *transcript_id) 
{
	char **tmp = NULL;

        if (*del == NULL) {
                if ((*del = calloc(1, sizeof(char *))) == NULL) {
                        printf("memory allocation error in add_deleted_transcript()\n");
                        exit(0);
                }
		if (((*del)[0] = calloc(strlen(transcript_id) + 1, sizeof(char))) == NULL) {
                        printf("memory allocation error in add_deleted_transcript()\n");
                        exit(0);
                }
		strcpy((*del)[0], transcript_id);
        }
        else {
                if ((tmp = (char **)realloc(*del, ((old + 1) * sizeof(char *)))) == NULL) {
                        printf("memory allocation error in add_deleted_transcript - realloc\n");
                        exit(0);
                }
                *del = tmp;
		if (((*del)[old] = calloc(strlen(transcript_id) + 1, sizeof(char))) == NULL) {
                        printf("memory allocation error in add_deleted_transcript\n");
                        exit(0);
                }
		strcpy((*del)[old], transcript_id);
        }
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void alter_description(char ***header, int cnt, char *vcf)
{
	int i,j;
	char tmp[512];

	tmp[0] = '\0';
	for (i=0; i<cnt; ++i) {
		if (strstr((*header)[i], "description") != NULL) {
			while ((*header)[i][(strlen((*header)[i] +1))] == '\r' || (*header)[i][(strlen((*header)[i] +1))] == '\n')
				(*header)[i][(strlen((*header)[i] +1))] = '\0';
			sprintf(tmp, "%s, liftover from %s\n", (*header)[i], vcf);
			if (tmp[0] != '\0') {
				free((*header)[i]);
				if (((*header)[i]=calloc((strlen(tmp) + 1), sizeof(char)))==NULL) {
                                	printf("memory allocation error in alter_description\n");
                                	exit(0);
				}
				strcpy((*header)[i], tmp);
			}
		}
	}
}
/************************************************************************************************/

/************************************************************************************************/
void parse_transcript_id(gtf_info *gtf)
{
	char *ptr = NULL, tmp[50];

	tmp[0] = '\0';
	if ((ptr = strstr((*gtf).group, "transcript_id")) != NULL) {
		sscanf(ptr + 13, "%s", (*gtf).transcript_id); 
	}
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void write_alt_gtf(FILE *f, gtf_info *gtf_lines, int gtf_cnt, contig_info contig)
{
	/*Delete whole transcript if start or end pos has been deleted:
		-parse gtf_lines to find transcript_id
		-scroll through gtf_lines, make list of deleted transcript_ids for each contig
		-scroll through gtf_lines, flag 'dont_print' if include deleted transcript_ids.
	For each gtf line, if !dont_print, look at st and en.
	if has transcript_id, check if st -> en includes any mut.pos. If so, add mut.mutation to gtf_lines group 
	Need to add the correct offset to each from contig.mut; find correct mut by eg. st > mut.pos[0] && st < mut.pos[+1]; add offset.
	*/
	int i = 0, j, k, l, flag = 0, index, del_transcript_cnt = 0;
	char **deleted_transcript = NULL;

	for (i=0; i<gtf_cnt; ++i) {
		parse_transcript_id(&gtf_lines[i]);
	}
	for (j = 0; j < gtf_cnt; ++j) {
		for (k = 0; k <contig.del_cnt; ++k) {
			if (gtf_lines[j].transcript_id && (gtf_lines[j].st == contig.deleted_pos[k] || gtf_lines[j].en == contig.deleted_pos[k])) {
				for (l = 0; l < del_transcript_cnt; ++l) {
					if (strcmp(deleted_transcript[l], gtf_lines[j].transcript_id) == 0)
					flag = 1; 
				}	
				if (!flag) {					
					add_deleted_transcript(&deleted_transcript, del_transcript_cnt, gtf_lines[j].transcript_id);
					del_transcript_cnt++;
				}
				flag = 0; 
			}
		}
	}
	for (j = 0; j < gtf_cnt; ++j) 
		for (k = 0; k < del_transcript_cnt; ++k)
			if (strcmp(gtf_lines[j].transcript_id, deleted_transcript[k]) == 0) 	
				gtf_lines[j].dont_print = 1;
	printf("deleted transcripts for contig %s:\n",contig.id);
	for (j=0; j<del_transcript_cnt; ++j) {
		printf("%d: %s\n",j,deleted_transcript[j]);
	} 
	for (j = 0; j < gtf_cnt; ++j) 
		if (gtf_lines[j].dont_print)
			printf("gtf omits line %d because %s was deleted\n", j, gtf_lines[j].transcript_id);	
//For each gtf line, if !dont_print, match seqname to contig.id
	flag = 0;	
	for (j = 0; j < gtf_cnt; ++j) { 
		if (!gtf_lines[j].dont_print) {
//if has transcript_id, check if st -> en includes any mut.pos. If so, add mut.mutation to gtf_lines group 
			if (gtf_lines[j].transcript_id[0] != '\0') {
				for (k = 0; k<contig.mut_cnt; ++k) {
					if ((gtf_lines[j].st -1 < contig.mut[k].pos) && (gtf_lines[j].en +1 > contig.mut[k].pos)) {
						contig.mut[k].flag = 1;
						flag = 1;													
					}
				}
				if (flag) {
					alter_group(gtf_lines[j].group, contig.mut, contig.mut_cnt);
					reset_flags(contig.mut, contig.mut_cnt);
					flag = 0;
				}
			}
//add the correct offset to each st and en from contig.mut
			index = 0;
			while ((index < contig.mut_cnt && contig.mut[index].pos < gtf_lines[j].st)) {
				index++;
			}
			if (index != 0) {	
				gtf_lines[j].st += contig.mut[index-1].offset;
			}
			index = 0;
			while ((index < contig.mut_cnt && contig.mut[index].pos < gtf_lines[j].en)) {
				index++;
			}
			if (index != 0) {
				gtf_lines[j].en += contig.mut[index-1].offset;
			}
		}
	}
      	for (i=0; i<gtf_cnt; i++) {
        	if (!gtf_lines[i].dont_print)
		        fprintf(f, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\t%s\n", gtf_lines[i].seqname, gtf_lines[i].source, gtf_lines[i].feature, gtf_lines[i].st, gtf_lines[i].en, gtf_lines[i].score, gtf_lines[i].strand, gtf_lines[i].frame, gtf_lines[i].group);
	}
	for (j=0; j<del_transcript_cnt; ++j) {
		if (deleted_transcript[j] != NULL)
			free(deleted_transcript[j]);
	}
	if (deleted_transcript != NULL)
		free(deleted_transcript);	
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void write_norm_gtf(FILE *f, gtf_info *gtf_lines, int gtf_cnt)
{
	int i;

      	for (i=0; i<gtf_cnt; i++) {
		fprintf(f, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\t%s\n", gtf_lines[i].seqname, gtf_lines[i].source, gtf_lines[i].feature, gtf_lines[i].st, gtf_lines[i].en, gtf_lines[i].score, gtf_lines[i].strand, gtf_lines[i].frame, gtf_lines[i].group);
	}
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void alter_group(char *group, mut_info *mut, int cnt) 
{
	int i;

	strcat(group, " mutations \"");
	for (i=0; i<cnt; ++i) {
		if (mut[i].flag) {
			strcat(group, mut[i].mutation);
			strcat(group, ".");
		}
	}
	if ((group)[(strlen(group) -1)] == '.')
		(group)[(strlen(group) -1)] = '\0';
	strcat(group, "\";");
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void reset_flags(mut_info *m, int cnt)
{
	int i;
	for (i=0; i<cnt; ++i)
		m[i].flag = 0;
	return;
}
/************************************************************************************************/

/************************************************************************************************/
void read_and_write_gtf(FILE *g, FILE *h, contig_info *contigs, int contig_cnt, gtf_metadata *gtf, int gtf_cnt)
{
	gtf_info *tmp = NULL;
        int i = 0, j = 0, k=0, match = -1;
        char line[MAXSTR], *pnt = NULL, *ptr = NULL;

	for (i = 0; i < gtf_cnt; ++i) {
		match = -1;
		tmp = NULL;
		if ((tmp=calloc(gtf[i].cnt, sizeof(gtf_info)))==NULL) {
               		printf("memory allocation error in read_and_write_gtf\n");
               		exit(0);
		}
		if (j != 0) {
                        sscanf(line, "%[^\t] %[^\t] %[^\t] %d %d %s %c %c %[^\r\n]", tmp[0].seqname, tmp[0].source, tmp[0].feature, &tmp[0].st, &tmp[0].en, tmp[0].score, &tmp[0].strand, &tmp[0].frame, tmp[0].group);
			j = 1;
		}
		
		while(fgets(line,MAXSTR - 1,g) != NULL && j < gtf[i].cnt) { 
			if (strncmp(line, gtf[i].seqname, strlen(gtf[i].seqname)) == 0) {
                        	sscanf(line, "%[^\t] %[^\t] %[^\t] %d %d %s %c %c %[^\r\n]", tmp[j].seqname, tmp[j].source, tmp[j].feature, &tmp[j].st, &tmp[j].en, tmp[j].score, &tmp[j].strand, &tmp[j].frame, tmp[j].group);
                        	j++;
                	}
		}
		for (k = 0; k<contig_cnt; ++k) {
			if (strcmp(contigs[k].id, gtf[i].seqname) == 0)
				match = k; 
		}
		if (j != gtf[i].cnt) {
			printf("Problem reading gtf file for contig = %s, j = %d, gtf[i].cnt = %d\n", gtf[i].seqname, j, gtf[i].cnt);
			exit(0);
		}
		if (match <0 || contigs[match].mut_cnt == 0) {//need to simply write these gtf lines as they stand!
			printf("No match for this contig from gtf: %s. These will be written as they are in the *_alt.gtf\n", gtf[i].seqname);
			write_norm_gtf(h, tmp, gtf[i].cnt);			
		}
		else
			write_alt_gtf(h, tmp, gtf[i].cnt, contigs[match]);
		if (tmp != NULL)
			free(tmp);
        }
	return;
}
/************************************************************************************************/
