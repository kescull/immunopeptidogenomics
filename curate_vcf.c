#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#define MAXSTR 51200

typedef struct vcf_info vcf_info;

struct vcf_info {
        char chrom[100];
	int pos;
	char id[20];
	char ref[512];
	char alt[512];
	char qual[10];
	char filter[512];
	char info[5120];
	int is_complex;
        int dont_print;
	int matches;
};

int count_vcf(FILE *, int *);
vcf_info *read_vcf(FILE *, char ***, int, int);
void print_vcf_line(FILE *, vcf_info);
void remove_duplicates(vcf_info *, int);
void check_indels(vcf_info *, int, int);
void remove_complex(vcf_info *, int);
int is_complex(char *); 

/* curate_vcf: curates vcf file to provide unambiguous input.vcf containing simple substitutions and indels for FastaAlternateReferenceMaker.
Input: .vcf file. Output: curated .vcf file named (input)_curated.vcf. 
Option: -d	use this command-line option to 'deprioritise deletions' i.e. remove deletion mutations which would 
		remove downstream mutation sites
Save stdout as log file - provides detailed information on curated mutations */
/***************************************************************************************************************************/
void main(int argc,char **argv) 
{
	FILE *f, *g;
	static char **header = NULL; 
	char *pnt = NULL, stem[100], output[100];
	vcf_info *vcf_lines = NULL;
	int header_cnt = 0, lines_cnt = 0, i, deprioritise_del = 0, c;

	while ((c = getopt (argc, argv, "d")) != -1) {
		switch (c) {
			case 'd':
				deprioritise_del = 1;
				break;
			case '?':
				printf("Unknown option '-%c'.\n", optopt);
				exit(0);
			default:
				exit(0);
		}
	}
	printf("deprioritise_del = %d\n", deprioritise_del);
	if (optind != (argc - 1)) {
                printf("error: specify .vcf file to curate\n");
                exit(0);
        }

	if ((f = fopen(argv[optind],"r")) == NULL) {
                printf("Can't open file %s\n",argv[optind]);
                exit(0);
        }
	strcpy(stem, argv[optind]);
        if ((pnt = strstr(stem, ".vcf")) != NULL)
                *pnt = '\0';
        else {
                printf("error: specify .vcf file - use .vcf extension\n");
                exit(0);
        }
	if (deprioritise_del)
		sprintf(output, "%s_unmasked.vcf", stem);
	else
		sprintf(output, "%s_indel.vcf", stem);

	if ((g = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n",output);
                exit(0);
        }
	
	lines_cnt = count_vcf(f, &header_cnt);
	printf("Found %d header lines and %d lines of mutation information\n", header_cnt, lines_cnt);
	
	vcf_lines = read_vcf(f, &header, lines_cnt, header_cnt);

	remove_complex(vcf_lines, lines_cnt);
	remove_duplicates(vcf_lines, lines_cnt);
	check_indels(vcf_lines, lines_cnt, deprioritise_del);

	for (i=0; i<header_cnt; ++i) {
		if (header[i] != NULL)
			fprintf(g, "%s", header[i]);
	}
	for (i=0; i<lines_cnt; i++) {
		if (!vcf_lines[i].dont_print) {
			print_vcf_line(g, vcf_lines[i]);
		}
	}

	fclose(f);
	fclose(g);
	for (i=0; i<header_cnt; ++i) {
		if (header[i] != NULL)
			free(header[i]);
	}
	if (header != NULL)
		free(header);
	if (vcf_lines != NULL)
		free(vcf_lines);
	printf("Finished\n");
	exit(0);
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
int count_vcf(FILE *f, int *header_cnt)
{
	int cnt = 0;
	char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#')
                        *header_cnt+=1;
		else
			cnt++;
        }
        rewind(f);

        return(cnt);
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
vcf_info *read_vcf(FILE *f, char ***header, int lines_cnt, int header_cnt)
{
	vcf_info *tmp = NULL, *big_var = NULL;
	int i = 0, j = 0, var_len = 0;
	char line[MAXSTR], chrom[150],id[25],ref[515],alt[515],qual[15],filter[515],info[5130];	

	if ((*header=(char **)calloc(header_cnt, sizeof(char *)))==NULL) {
                printf("memory allocation error in read_vcf\n");
                exit(0);
        }
	if ((tmp=(vcf_info *)calloc(lines_cnt, sizeof(vcf_info)))==NULL) {
                printf("memory allocation error in read_vcf\n");
                exit(0);
        }
	
	while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '#' && i < header_cnt) {
                        if (((*header)[i]=(char *)calloc((strlen(line) + 1), sizeof(char)))==NULL) {
                                printf("memory allocation error in read_vcf\n");
                                exit(0);
                        }
                        strcpy((*header)[i], line);
                        i++;
                }
		else if (j < lines_cnt) {
			chrom[0] = id[0] = ref[0] = alt[0] = qual[0] = filter[0] = info[0] = '\0';
			sscanf(line, "%[^\t] %d %[^\t] %[^\t] %[^\t] %s %s %s", chrom, &tmp[j].pos, id, ref, alt, qual, filter, info);
			if (strlen(chrom) > 99 || strlen(id) > 19 || strlen(ref) > 511 || strlen(alt) > 511 || 
					strlen(qual) > 9 || strlen(filter) > 511 || strlen(info) > 5119) {
				printf("Content of one of the vcf fields is too big for the array specified in vcf_info structure. Please check your vcf, increase the size of the appropriate array and try again.\n");
				exit(0);
			}
			strcpy(tmp[j].chrom,chrom);
			strcpy(tmp[j].id,id);
			strcpy(tmp[j].ref,ref);
			strcpy(tmp[j].alt,alt);
			strcpy(tmp[j].qual,qual);
			strcpy(tmp[j].filter,filter);
			strcpy(tmp[j].info,info);

			if (strlen(tmp[j].ref) > var_len) {
				var_len = strlen(tmp[j].ref);
				big_var = &tmp[j];
			}
			if (strlen(tmp[j].alt) > var_len) {
				var_len = strlen(tmp[j].alt);
				big_var = &tmp[j];
			}
			j++;
		}
		else {
			printf("Not enough lines in vcf_lines\n");
			exit(0);
		}
	}
	printf("biggest variant is %d:\n", var_len);
	printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", big_var->chrom, big_var->pos, big_var->id, big_var->ref, big_var->alt, big_var->qual, big_var->filter, big_var->info);
	
	return tmp;
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
void print_vcf_line(FILE *f, vcf_info v)
{
	fprintf(f, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", v.chrom, v.pos, v.id, v.ref, v.alt, v.qual, v.filter, v.info);
	return;
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
void remove_duplicates(vcf_info *v, int cnt)
{
	int i, j = 1, k, l, highest = 0, flag = 0;

	printf("The following sets of mutations were found at the same genomic coordinates (some may be exact duplicates).\n");
	printf("Only one of each set has been retained in the curated .vcf\n\n");

	for (i=0; i<cnt; ++i) {
		if (!v[i].is_complex) {
			while(i + j < cnt && (strcmp(v[i].chrom, v[i+j].chrom) == 0) && (v[i].pos == v[i+j].pos)) { //count how many in a row with matching chrom and pos (=> j)
				++j;
			}
			if (j > 1) {
				//printf("Found multiple mutations at this site:\n");
				for (k=0; k<j-1; k++) {
					for (l=1; l<j-k; l++) { //for each in the set, check how many exactly match it
						if (!v[i+k].is_complex && !v[i+k+l].is_complex)
							if (strcmp(v[i+k].ref,v[i+k+l].ref) == 0 && strcmp(v[i+k].alt,v[i+k+l].alt) == 0) {
								v[i+k].matches++;
								v[i+k+l].matches++;
							}	
					}
				}
				for (k=i; k<i+j; k++) { //print a log of multiples and their matches, and record highest number of matches for the set
					if (!v[k].is_complex) {
						printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\tmatches %d\n", v[k].chrom, v[k].pos, v[k].id, v[k].ref, v[k].alt, v[k].qual, v[k].filter, v[k].info, v[k].matches);
						if (v[k].matches > highest)
							highest = v[k].matches;
						v[k].dont_print = 1;
					}
				}
				for (k=i; k<i+j; k++) {
//need new flag! dont_print is already flagged for 'complex' now!
					if (!v[k].is_complex)
						if (v[k].matches < highest) { //print to log any mutations that won't be represented in output.vcf
							printf("This mutation will not be in the curated .vcf because there is a more abundant mutation at this position:\n");
							printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", v[k].chrom, v[k].pos, v[k].id, v[k].ref, v[k].alt, v[k].qual, v[k].filter, v[k].info);
							//fflush(stdout);
						}
					if (!flag) { //ensure only the first of those with the highest number of matches is printed in output.vcf
						if (v[k].matches == highest) {
							v[k].dont_print = 0;
							flag=1;
						}
					}
				}
				printf("\n");
			}			
			i+=j-1;
			j=1;
			highest=0;
			flag=0;
		}
	}
	return;
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
void check_indels(vcf_info *v, int cnt, int deprioritise_del) 
{
	int i, j,del_cnt=0, removed_cnt = 0, offset = 0, flag = 0;

	for (i = 0; i < cnt; ++i) {
		if (!v[i].dont_print) {
			offset = strlen(v[i].ref) - strlen(v[i].alt);
			if (offset > 0) { //deletion mutation
				j=1;
				del_cnt++;
				if (deprioritise_del) {
					while ((strcmp(v[i].chrom, v[i+j].chrom) == 0) && (v[i+j].pos < (v[i].pos + offset + 1))) {
						if (!v[i+j].dont_print) {
							flag = 1;
						}
						j++;
					}	
					if (flag) {
						v[i].dont_print = 1;
						printf("This deletion mutation will not be in the curated .vcf because including it would remove downstream mutation sites:\n");
						printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", v[i].chrom, v[i].pos, v[i].id, v[i].ref, v[i].alt, v[i].qual, v[i].filter, v[i].info);
						removed_cnt++;
					}
					flag = 0;		
				}
				else {
					while ((strcmp(v[i].chrom, v[i+j].chrom) == 0) && (v[i+j].pos < (v[i].pos + offset + 1))) {
						if (!v[i+j].dont_print) {
							v[i+j].dont_print = 1;
							printf("This mutation will not be in the curated .vcf because a deletion of %d bases at chrom %s, pos %d has removed this position:\n", offset, v[i].chrom, v[i].pos);
							printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", v[i+j].chrom, v[i+j].pos, v[i+j].id, v[i+j].ref, v[i+j].alt, v[i+j].qual, v[i+j].filter, v[i+j].info);
							removed_cnt++;
						}	
						j++;
					}
				}			
			}
		}
	}
	printf("Found %d deletions, resulting in %d lines removed\n", del_cnt, removed_cnt);
	if (!removed_cnt)
		printf("No lines were deleted - the output vcf is identical to the original\n");

	return;
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
void remove_complex(vcf_info *v, int cnt)
{
	int i, j, flag = 0;
	
	printf("The following mutations have been removed because they are too complex for FastaAlternateReferenceMaker:\n");

	for (i = 0; i < cnt; ++i) {
		if (is_complex(v[i].ref) || is_complex(v[i].alt)) {
			v[i].is_complex = 1;
			v[i].dont_print = 1;
			printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", v[i].chrom, v[i].pos, v[i].id, v[i].ref, v[i].alt, v[i].qual, v[i].filter, v[i].info);
		} 
	}
	printf("\n");
	return;
}
/***************************************************************************************************************************/

/***************************************************************************************************************************/
int is_complex(char *s) 
{
	char *simple = "GATC";
	int i = 0;
	for (i=0; i<strlen(s); ++i) {
		if (strchr(simple, s[i]) == NULL)
			return 1;
	}
	return 0;
}
/***************************************************************************************************************************/
