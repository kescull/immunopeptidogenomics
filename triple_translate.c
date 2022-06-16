#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#define MAXSTR 51200
#define FASTALEN 70
#define MAXSEQLEN 350000
#define MINLEN 8

typedef enum {
	INDEL,
	SUB,
	NO_CHANGE,
} Mut;

typedef enum {
	PLUS,
	MINUS,
	UNKNOWN,
} Strand;
typedef enum {
	DEFAULT,
	CUFFCOMPARE,
} mode_type;

typedef struct mut_info mut_info;
typedef struct seq_info seq_info;
typedef struct exon exon;
typedef struct aa_mut_info aa_mut_info;
typedef struct gene_info gene_info;

struct gene_info {
	char transcript[512];
	char gene[512];
	char class_code;
};

struct exon {
	int st;
	int en;
	int len;
	int transcript_st;
	int transcript_en;
};
struct aa_mut_info {
	Mut type;
	char desc[5210];
	int pos;
};
struct mut_info {
	int coordinate;
	int seq_loc;
	char from[512];
	char to[512];
	aa_mut_info aa_muts[3];
};
struct seq_info {
	int mut_cnt;
	mut_info *muts;
	Strand strand;
	char id[50];
	char loc[50];
	char chr[50];
	char gene[512];
	char class_code;
};

mode_type mode = DEFAULT; 
gene_info *genes = NULL;
int gene_cnt = 0;

void translate_fasta(FILE *, FILE *);
void print_fasta_seq(FILE *, char *, char *);
char **make_translation(char *, seq_info *info);
seq_info parse_header(char *);
char *move_ptr(char *,char);
int parse_mutations(mut_info **, char *);
int parse_exons(exon **, char *, char *, int);
void get_mut_loc(exon *, int, mut_info *, int, int);
void format_translation(FILE *, char **, seq_info);
int count_genes(FILE *);
gene_info *parse_gene_info(FILE *, int);
void print_usage();
/* triple_translate: take transcriptome nucleotide fasta file and perform 3 frame translation.
Keep any sequences > 7 aa in length. Label by transcript name, frame, any mutations
Output: protein fasta file for searching with PEAKS etc
 */
/*************************************************************************************************/
void main(int argc,char **argv)
{
	FILE *f = NULL, *g = NULL, *h = NULL;
	char stem[100], output[100], *pnt = NULL, c, *cvalue = NULL;
	int i;

	mode = DEFAULT;
	while ((c = getopt (argc, argv, "c:h")) != -1) {
                switch (c) {
                        case 'c':
				printf("Cuffcompare options triggered\n");
				mode = CUFFCOMPARE;
                                cvalue = optarg;
				if (strstr(optarg, ".tracking") == NULL) {
          				printf ("Specify the cuffcompare .tracking file following option -c\n\n");
					print_usage();
					exit(0);
				}
                                break;
			case 'h':
				print_usage();
				exit(0);
                        default: /*'?'*/
				if (optopt == 'c')
          				printf ("Specify the cuffcompare .tracking file following option -c\n\n");
				else
                                	printf("Unknown option '-%c'.\n\n", optopt);
				print_usage();
                                exit(0);
                }
        }	
	
	if (optind != (argc - 1)) {
		printf("Specify transcriptome .fa/.fasta file for translation\n\n");
		print_usage();
                exit(0);
        }

	if ((f = fopen(argv[optind],"r")) == NULL) {
                printf("Can't open file %s\n",argv[optind]);
                exit(0);
        }
	if (mode == CUFFCOMPARE) {
		if ((h = fopen(cvalue,"r")) == NULL) {
        	        printf("Can't open file %s\n",cvalue);
        	        exit(0);
        	}
		gene_cnt = count_genes(h);
		genes = parse_gene_info(h, gene_cnt);
		printf("Transcripts with associated known genes: %d\n", gene_cnt);
	}
	
	strcpy(stem, argv[optind]);
	if ((pnt = strstr(stem, ".fa")) != NULL) {
		*pnt = '\0';
	}
	sprintf(output, "%s_3translate.fasta", stem);
	if ((g = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n",output);
                exit(0);
        }

	translate_fasta(f, g);

	if (mode == CUFFCOMPARE) {
		if (genes != NULL)
			free(genes);
		fclose(h);
	}
		
	fclose(f);
	fclose(g);

	exit(0);
}
/*************************************************************************************************/

/*************************************************************************************************/
void print_usage()
{
	printf("Required input:\n\tspecify one transcriptome sequences .fa/.fasta file for translation\n");
	printf("Optional:\n\t-c\tCuffcompare output tracking file\n");
	printf("\t-h\tprint usage help\n");
	return;
}
/*************************************************************************************************/

/*************************************************************************************************/
void translate_fasta(FILE *f, FILE *g)
{
	int transcript_cnt = 0, len = 0, first = 1, i=0, frame, unknown_strand_cnt = 0;
	char line[MAXSTR], header[MAXSTR], seq[MAXSEQLEN], **aa = NULL;
	seq_info info;

	header[0] = '\0';
	seq[0] = '\0';
	
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (line[0] == '>') {
			if (first) {
				strcpy(header, line);
				info = parse_header(header);
				first = 0;
			}
			else {
				if (info.strand == UNKNOWN) {
					unknown_strand_cnt++;
					printf("Transcript with unknown strand not printed:\n\t%s\n", header);
				} 
				else {
				//deal with seq
					aa = make_translation(seq, &info);
					format_translation(g, aa, info);
					for (i = 0; i<3; i++) {
						if (aa[i] != NULL)		
							free(aa[i]);
					}
					if (aa != NULL)
						free(aa);
				}
				transcript_cnt++;
//empty the variables
				seq[0] = '\0';
				len = 0;
				header[0] = '\0';
				if (info.mut_cnt >0) {
					if(info.muts != NULL)
						free(info.muts);
				}
//Start new transcript
				strcpy(header,line);
				info = parse_header(header);
			}
		}
		else {
			while (line[strlen(line) - 1] == '\n' || line[strlen(line) - 1] == '\r')
				line[strlen(line) - 1] = '\0';
			len += strlen(line);
			if (len < MAXSEQLEN)
				strcat(seq, line);
			else {
				printf("sequence exceeds max length of %d - increase MAXSEQLEN\n", MAXSEQLEN);
				exit(0);
			}
		}
	}
//last
	//deal with seq	
	if (info.strand == UNKNOWN) {
		unknown_strand_cnt ++;
		printf("Transcript with unknown strand not printed:\n\t%s\n", header);
	} else { 
		aa = make_translation(seq, &info);
		for (i = 0; i<3; i++) {
			if (aa[i] != NULL)
				free(aa[i]);
		}
		if (aa != NULL)
			free(aa);
	}
	transcript_cnt++;
	printf("Transcripts: %d; Unknown strand: %d\n", transcript_cnt, unknown_strand_cnt);
	if (info.mut_cnt >0) 
		if(info.muts != NULL)
			free(info.muts);
	return;
}
/*************************************************************************************************/

/*************************************************************************************************/
void print_fasta_seq(FILE *f, char *seq, char *header)
{
	int len = strlen(seq), i = 0, end=0;
	char copy[FASTALEN + 2];

	fprintf(f, "%s", header);
	while (i < len) {
		copy[0] = '\0';
		strncpy(copy, &seq[i], FASTALEN);
		copy[FASTALEN] = '\0';
		end = strlen(copy);	
		copy[end] = '\n';
		copy[end+1] = '\0';
		fprintf(f, "%s", copy);
		i += FASTALEN;
	}
}
/*************************************************************************************************/

/*************************************************************************************************/
char **make_translation(char *na, seq_info *info)
{
	int x,y,len, i, j,k, pos, frame, index = 0, shift = 0, codon = 0;
	char c[5][5][5] = {{"FFLLX","SSSSX","YYZZX","CCZWX","XXXXX"},  //N.B.: Z is STOP
                   	{"LLLLX","PPPPX","HHQQX","RRRRX","XXXXX"},
                   	{"IIIMX","TTTTX","NNKKX","SSRRX","XXXXX"},
                   	{"VVVVX","AAAAX","DDEEX","GGGGX","XXXXX"},
                   	{"XXXXX","XXXXX","XXXXX","XXXXX","XXXXX"}};
	int coded[] = {2,4,1,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,0,0,4,4,4,4,4};
	int na_len = strlen(na);
	char **aa = NULL, old[4], new[4], aa_new[2], aa_old[2];
	
	if ((aa = (char **)calloc(3,sizeof(char *))) == NULL) {
                printf("Memory allocation error in make_translation()\n");
                exit(0);
        }
        for (x = 0;x < 3;x ++) {
                if((aa[x] = (char *)calloc(((na_len) / 3) + 3,sizeof(char))) == NULL) {
                        printf("Memory allocation error in search_na_sequence()\n");
                        exit(0);
                }
        }
	for (i = 0; i< info->mut_cnt; i ++) {
		if (info->muts[i].seq_loc > -1) {
			if ((strlen(info->muts[i].to) > 1) || (strlen(info->muts[i].from) > 1)) {
				for (frame=0; frame<3; frame++) {
					info->muts[i].aa_muts[frame].type = INDEL;
					info->muts[i].aa_muts[frame].desc[0] = '\0'; 
					sprintf(info->muts[i].aa_muts[frame].desc, "INDEL:(%s->%s@%s:%d)", info->muts[i].from, info->muts[i].to, info->chr, info->muts[i].coordinate);
					info->muts[i].aa_muts[frame].pos = (info->muts[i].seq_loc - frame)/3;
				}
			}
			else {
				new[0] = '\0';
				old[0] = '\0';
				index = info->muts[i].seq_loc -1;
		//for each frame, work out which codon type (x in start, middle, or end) and codon position.
		//work out codon translations, and see if change. if so, save in aa_muts for that frame - type SUB 
				for (frame = 0; frame < 3; frame++) {
					pos = 0;
					shift = (index - frame)%3;
					codon = (index - frame)/3;
					for (k=0;k<3;++k) {
						pos = index - shift + k;
						if (pos > -1 && pos < na_len) {
							new[k] = na[pos];
							if ((k - shift) == 0) {
								if (info->strand == PLUS) {
									old[k] = *info->muts[i].from;
								}	
								else if (info->strand == MINUS) {
									switch (*info->muts[i].from) {
										case 'A':
											old[k] = 'T';
											break;
										case 'T':
											old[k] = 'A';
											break;
										case 'C':
											old[k] = 'G';
											break;
										case 'G':
											old[k] = 'C';
											break;
										default:
											printf("Found a nucleotide that is not GATC %s\n", info->muts[i].from);
											exit(0);	
									}
								}
								else {
									printf("Unknown strand in make_translation()\n");
									exit(0);
								}
							} else {
								old[k] = na[pos];
							}
						}
					}	
					old[k] = '\0';
					new[k] = '\0';
					aa_old[0] = '\0'; 
					aa_new[0] = '\0'; 
					if (strlen(old) == 3 && strlen(new) == 3) {
                        			aa_old[0] = c[coded[old[0] - 'A']][coded[old[1] - 'A']][coded[old[2] - 'A']];
                        			aa_new[0] = c[coded[new[0] - 'A']][coded[new[1] - 'A']][coded[new[2] - 'A']];
						aa_old[1] = '\0'; 
						aa_new[1] = '\0'; 
						if (strcmp(aa_new, aa_old) != 0) {
							info->muts[i].aa_muts[frame].type = SUB;
							sprintf(info->muts[i].aa_muts[frame].desc, "SUB:%s->%s", (aa_old[0] == 'Z') ? "STOP" : aa_old, (aa_new[0] == 'Z') ? "STOP" : aa_new);
							info->muts[i].aa_muts[frame].pos = codon;
						} else {
							info->muts[i].aa_muts[frame].type = NO_CHANGE;
						}
					} else {
						info->muts[i].aa_muts[frame].type = NO_CHANGE;
					}
				}
			}
		}
	}
	for(x = 0;x < 3;x ++) {
        	len = 0;
        	for(y = x;y < na_len;y += 3) {
                	if(na[y] != '\0' && na[y + 1] != '\0' && na[y + 2] != '\0') {
                        	aa[x][len] = c[coded[na[y] - 'A']][coded[na[y + 1] - 'A']][coded[na[y + 2] - 'A']];
                        	len ++;
                        	aa[x][len] = '\0';
                        }
                }
        }

return(aa);
}
/*************************************************************************************************/

/*************************************************************************************************/
seq_info parse_header(char *header)
{
	seq_info tmp;
	exon *exons = NULL;
	char *ptr = NULL, *pnt = NULL, chunk[strlen(header)], chunk2[strlen(header)], c = '\0';
	int i, exon_cnt, flag = 0;
	
	chunk[0] = '\0';
	chunk2[0] = '\0';
	tmp.chr[0] = '\0';

//get id and loc
	sscanf(header+1, "%s", tmp.id);
	ptr = strstr(header, "loc:");
	sscanf(ptr, "%s", tmp.loc);	
	
//get strand
	c = tmp.loc[strlen(tmp.loc)-1];
	if (c == '+')
		tmp.strand = PLUS;
	else if (c == '-')
                tmp.strand = MINUS;
	else {
		tmp.strand = UNKNOWN;
	}
//get chr
	ptr = tmp.loc + 4;
	pnt = move_ptr(ptr, '|');
	strcpy(tmp.chr, ptr);
	*pnt = '|';

//get mutations
	tmp.mut_cnt = 0;
	if ((ptr = strstr(header, ";mutations=")) != NULL) {
		strcpy(chunk, ptr+11);
		pnt = move_ptr(ptr, ';');
		tmp.mut_cnt = parse_mutations(&tmp.muts, chunk);
		if (pnt != NULL)
			*pnt = ';';
	}
	chunk[0] = '\0';
	chunk2[0] = '\0';
	if (tmp.mut_cnt > 0 && tmp.strand != UNKNOWN) {
		if ((ptr = strstr(header, "exons:")) != NULL) {
			ptr += 6;
			pnt = move_ptr(ptr, ' ');
			strcpy(chunk, ptr);
			*pnt = ' ';
			ptr = pnt + 6;
			pnt = move_ptr(ptr, ' ');
			strcpy(chunk2, ptr);
			*pnt = ' ';
			exon_cnt = parse_exons(&exons, chunk, chunk2, tmp.strand);
			if (exon_cnt > 0) {
				get_mut_loc(exons, exon_cnt, tmp.muts, tmp.mut_cnt, tmp.strand);

			} else {
				printf("parse_exons failed...\n");
				exit(0);
			}	
		}  
	}
	if (exons != NULL)
		free(exons);
//if CUFFCOMPARE, match gene if there is one
	tmp.gene[0] = '\0';
	if (mode == CUFFCOMPARE) {
		for (i=0; i<gene_cnt; ++i) {
			if ((strcmp(tmp.id, genes[i].transcript)) == 0) {
				strcpy(tmp.gene, genes[i].gene);
				tmp.class_code = genes[i].class_code;
				break;
			}
		}
	}
	
	return tmp;
}
/*************************************************************************************************/

/*************************************************************************************************/
char *move_ptr(char *str,char delim)
{
char *tmp = NULL;

if((tmp = strchr(str,delim)) != NULL) *tmp = '\0';

return(tmp);
}
/*************************************************************************************************/

/*************************************************************************************************/
int parse_mutations(mut_info **muts, char *chunk)
{
	int cnt = 0, i=0;
	char *pnt = NULL, *ptr = NULL, copy[strlen(chunk) + 1];

	strcpy(copy, chunk);
	while (*chunk != '\n' && *chunk != ';') {
		if (*chunk++ == '@')
			cnt++;
	}
	if ((*muts=(mut_info *)calloc(cnt, sizeof(mut_info)))==NULL) {
                printf("memory allocation error in parse_mutations\n");
                exit(0);
        }
	ptr = copy;
	for (i = 0; i < cnt; i++) {
		pnt = move_ptr(ptr, '-');
		strcpy((*muts)[i].from, ptr);
		ptr = pnt +2;
		pnt = move_ptr(ptr, '@');
		strcpy((*muts)[i].to, ptr);
		ptr = pnt +1;
		pnt = (i == cnt - 1) ? move_ptr(ptr, '\n') : move_ptr(ptr, '.');
		(*muts)[i].coordinate = atoi(ptr);	
		if (i < cnt -1)
			ptr = pnt +1;
	}
	return cnt;
}
/*************************************************************************************************/

/*************************************************************************************************/
int parse_exons(exon **exons, char *chunk, char *chunk2, int strand)
{
	int cnt = 0, i=0;
	char *pnt = NULL, *ptr = NULL, copy[strlen(chunk) + 1];

	strcpy(copy, chunk);
	while (*chunk != '\0') {
		if (*chunk++ == '-')
			cnt++;
	}
	if ((*exons=(exon *)calloc(cnt, sizeof(exon)))==NULL) {
                printf("memory allocation error in parse_exons\n");
                exit(0);
        }
	ptr = copy;
//copy gene loci into exon structure array
	for (i = 0; i < cnt; i++) {
		pnt = move_ptr(ptr, '-');
		(*exons)[i].st = atoi(ptr);
		ptr = pnt +1;
		pnt = move_ptr(ptr, ',');
		(*exons)[i].en = atoi(ptr);	
		if (i < cnt -1)
			ptr = pnt +1;
	}		
	ptr = chunk2;
//copy transcript loci into exon structure array
	if (strand == PLUS) {
		for (i = 0; i < cnt; i++) {
			pnt = move_ptr(ptr, '-');
			(*exons)[i].transcript_st = atoi(ptr);
			ptr = pnt +1;
			pnt = move_ptr(ptr, ',');
			(*exons)[i].transcript_en = atoi(ptr);	
			if (i < cnt -1)
				ptr = pnt +1;
		}		
	} else if (strand == MINUS) {
		for (i = cnt-1; i > -1; i--) {
			pnt = move_ptr(ptr, '-');
			(*exons)[i].transcript_st = atoi(ptr);
			ptr = pnt +1;
			pnt = move_ptr(ptr, ',');
			(*exons)[i].transcript_en = atoi(ptr);	
			if (i > 0)
				ptr = pnt +1;
		}		

	} else {
		printf("Strand is neither PLUS nor MINUS in parse_exons()\n");
		exit(0);
	}
	for (i = 0; i < cnt; ++i) {
		(*exons)[i].len = (*exons)[i].en - (*exons)[i].st + 1;
		if ((*exons)[i].len != ((*exons)[i].transcript_en - (*exons)[i].transcript_st + 1)) {
			printf("Exon and seg length don't match!\n");
			exit(0);
		}
	}
	
	return cnt;
}
/*************************************************************************************************/
	
/*************************************************************************************************/
void get_mut_loc(exon *exons, int exon_cnt, mut_info *muts, int mut_cnt, int strand)
{
	int i, j;

	for (i=0; i<mut_cnt; ++i) {
		muts[i].seq_loc = -1;
		for (j=0; j<exon_cnt; ++j) {
			if ((muts[i].coordinate < (exons[j].en +1)) && (muts[i].coordinate > (exons[j].st -1))) {
				if (strand == PLUS) {
					muts[i].seq_loc = exons[j].transcript_en - exons[j].en + muts[i].coordinate;
				} else if (strand == MINUS) {
					muts[i].seq_loc = exons[j].transcript_en - muts[i].coordinate + exons[j].st;
				} else {
					printf("get_mut_loc was handed unstranded info\n");
					exit(0);
				}
			}
		}
	}
	return;
}
/*************************************************************************************************/

/*************************************************************************************************/
void format_translation(FILE *f, char **aa, seq_info info)
{
	int i, j, frame, pos = 0, tmp_pos = 0, seg_cnt = 0;
	char header[MAXSTR], seq[MAXSEQLEN], *pnt = NULL, *ptr = NULL, cat_desc[MAXSTR], tmp[MAXSTR], gene[520];

//Need to: divvy up each frame seq into chunks between 'STOP', check for mutation, make new headers, print.
	for (frame = 0; frame < 3; ++frame) {
		seq[0] = header[0] = cat_desc[0] = tmp[0] = gene[0] = '\0';	
		seg_cnt = 0;
		for (i=0, pos = 0; pos < strlen(aa[frame]); ++pos) {
			if (aa[frame][pos] == 'Z') {
				seq[i] = '\0';
				if (strlen(seq) > (MINLEN - 1)) {
					seg_cnt++;
					for (j=0; j<info.mut_cnt; j++) {
						if (info.muts[j].aa_muts[frame].type != NO_CHANGE) {
							if (info.muts[j].aa_muts[frame].pos < (pos+1) && info.muts[j].aa_muts[frame].pos > (pos-i)) {
								if ((tmp_pos = i - pos + info.muts[j].aa_muts[frame].pos + 1) == (i + 1))
									sprintf(tmp, " %s@END;",info.muts[j].aa_muts[frame].desc);
								else
									sprintf(tmp, " %s@%d;",info.muts[j].aa_muts[frame].desc, tmp_pos);
								strcat(cat_desc, tmp);
								tmp[0] = '\0';
							}
						}
					}
					if (mode == CUFFCOMPARE) {
						if (info.gene[0] != '\0') 
							sprintf(gene, " %s:%c", info.gene, info.class_code);
					}
					
					sprintf(header, ">%s_f%dp%d%s%s\n", info.id, frame + 1, seg_cnt, gene, cat_desc);
					print_fasta_seq(f, seq, header);
				}
//reset
				header[0] = cat_desc[0] = seq[0] = gene[0] = '\0';
				i = 0;
				
			}
			else {
				seq[i++] = aa[frame][pos];
			}
		}
//deal with last
		seq[i] = '\0';
		if (strlen(seq) > (MINLEN - 1)) {
			seg_cnt++;
			for (j=0; j<info.mut_cnt; j++) {
				if (info.muts[j].aa_muts[frame].type != NO_CHANGE) {
					if (info.muts[j].aa_muts[frame].pos < (pos+1) && info.muts[j].aa_muts[frame].pos > (pos-i)) {
						if ((tmp_pos = i - pos + info.muts[j].aa_muts[frame].pos + 1) == (i + 1))
							sprintf(tmp, " %s@END;",info.muts[j].aa_muts[frame].desc);
						else
							sprintf(tmp, " %s@%d;",info.muts[j].aa_muts[frame].desc, tmp_pos);
						strcat(cat_desc, tmp);
						tmp[0] = '\0';
					}
				}
			}
			if (mode == CUFFCOMPARE) {
				if (info.gene[0] != '\0') 
					sprintf(gene, " %s:%c", info.gene, info.class_code);
			}
			sprintf(header, ">%s_f%dp%d%s%s\n", info.id, frame + 1, seg_cnt, gene, cat_desc);
			print_fasta_seq(f, seq, header);
		}
	} 
	return;
}
/*************************************************************************************************/

/*************************************************************************************************/
int count_genes(FILE *f)
{
	int cnt = 0;
	char line[MAXSTR], tmp[MAXSTR];

	while(fgets(line,MAXSTR - 1,f) != NULL) {
		tmp[0] = '\0';
		sscanf(line, "%*s %*s %s", tmp);
		if (tmp[0] != '-')         
                	cnt++;
        }
        rewind(f);

        return(cnt);
}
/*************************************************************************************************/

/*************************************************************************************************/
gene_info *parse_gene_info(FILE *f, int cnt)
{
	int i=0;
	char line[MAXSTR], tmp_tran[512], tmp_gene[512], c;
	gene_info *tmp;

	if ((tmp = (gene_info *)calloc(cnt,sizeof(gene_info))) == NULL) {
                printf("Memory allocation error in parse_gene_info()\n");
                exit(0);
        }

	while(fgets(line,MAXSTR - 1,f) != NULL) {
		tmp_tran[0] = tmp_gene[0] = '\0';
		sscanf(line, "%s %*[^\t] %s %c", tmp_tran, tmp_gene, &c);
		if (tmp_gene[0] != '-') {
			strcpy(tmp[i].transcript, tmp_tran);
			strcpy(tmp[i].gene, tmp_gene);
			tmp[i].class_code = c;
			i++;	
		}
	}

	if (i != cnt) {
		printf("Counting mismatch in parse_gene_info()\n");
		exit(0);
	}
	return tmp;
}
/*************************************************************************************************/
