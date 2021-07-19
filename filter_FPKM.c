#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#define MAXSTR 51200

typedef enum {
        FPKM,
        CONF_HI,
} mode_type;

typedef struct {
	char gene[100];
	int found_in_rna;
} gene_evidence;

typedef struct {
	char transcript[100];
	char gene[100];
	int found_in_rna;
	int is_ref_transcript;
	int in_original_gtf;
} evidence;

void print_help();
int count_lines(FILE *);
evidence *get_good_transcripts(FILE *, double, int, mode_type);
char *move_ptr(char *, char);
char *skip_fields(char *, char, int);
int sort_gene_alpha(const void *, const void *); 
int is_good_transcript(char *, evidence **, int);

//filter_FPKM.c: finds list of transcripts with FPKM > 0 (or specified value)
//outputs new gtf filtered for transcripts that pass this test
//also outputs '*_table.csv' file with breakdown of info from tracking file,
//input and output gtfs. stdout contains various useful statistics. 
//
/************************************************************************/
void main(int argc, char **argv) 
{
	FILE *f, *g, *table_output;
	char c = '\0', gtf_fn[512], tracking_fn[512], threshold_string[100];
       	char output_fn[512], stem[512], *pnt = NULL, line[MAXSTR];
	double threshold = 0.0;
	int all_cnt = 0, good_cnt = 0, i,j, ref_cnt=0, good_ref_cnt=0, gene_cnt=0, good_gene_cnt = 0, max_gene_cnt = 0;
	int final_ref_cnt=0, final_good_ref_cnt=0, gene_nodup_cnt = 0, novel_t = 0, novel_good_t = 0, unfilt_gtf_t = 0, final_t = 0;
	mode_type mode = FPKM;
	evidence *evidence_table = NULL;
	gene_evidence *tracking_gene_table = NULL, *no_dup_tracking_gene_table = NULL;
	gene_evidence *gene_table = NULL, *no_dup_gene_table = NULL;

	gtf_fn[0] = tracking_fn[0] = threshold_string[0] = '\0';
	while ((c = getopt (argc, argv, "cg:t:f:h")) != -1) {
                switch (c) {
			case 'c':
				mode = CONF_HI;
                                printf("-a\tselect transcripts with any evidence (conf_hi > 0)\n");
                                break;
                        case 'g':
				strcpy(gtf_fn, optarg);
                                printf("-g\tgtf file %s\n", gtf_fn);
                                break;
			case 't':
				strcpy(tracking_fn, optarg);
                                printf("-t\ttracking file %s\n", tracking_fn);
                                break;
			case 'f':
				threshold = atof(optarg);
				printf("-f\tfilter threshold %.2lf\n", threshold);
				break;
			case 'h':
				print_help();
                                exit(0);
                        default: /*'?'*/
                                printf("Use options to specify files for analysis\n\n");
                                print_help();
                                exit(0);
                }
        }
	if (tracking_fn[0] == '\0' || gtf_fn[0] == '\0') {
		print_help();
		exit(0);
	}
	printf("\nIn the following statistics:\n\t'novel' is defined as having a cufflinks class code that is not = or c\n");
	printf("\tgenes are counted by gene_id\n");
	printf("\tcuffcompare removes some 'redundant' transcripts, so there are\n\t\tfewer transcripts in the unfiltered gtf than in the tracking file.\n");
	printf("\tFor calculating percentage coverage of reference transcripts/genes,\n\t\tthe unfiltered gtf is considered to contain 100 percent\n\n");

	//get good transcripts from tracking file (also calculate ref transcript and gene coverage)
	if ((f = fopen(tracking_fn,"r")) == NULL) {
        	printf("Can't open file %s\n",tracking_fn);
                exit(0);
        }
	all_cnt = count_lines(f);

	evidence_table = get_good_transcripts(f, threshold, all_cnt, mode);
	fclose(f);

	for (i = 0; i < all_cnt; ++i) {
		if(evidence_table[i].found_in_rna) {
			good_cnt++;
		} 
		if (evidence_table[i].is_ref_transcript) {
			ref_cnt++;
		} else {
			novel_t++;
			if (evidence_table[i].found_in_rna) {
				novel_good_t++;
			}
		}
		if (evidence_table[i].found_in_rna && evidence_table[i].is_ref_transcript) {
			good_ref_cnt++;
		}
	}
	sprintf(threshold_string, "%s > %.2lf",(mode == FPKM) ? "FPKM":"conf_hi",threshold);
	printf("Tracking file information:\n");
	printf("Total transcripts: %d\n",all_cnt);
	printf("Novel transcripts: %d (%.2lf%% of total transcripts)\n",novel_t, (double)novel_t/all_cnt*100);
	printf("Transcripts with %s: %d (%.2lf%% of total transcripts)\n",threshold_string, good_cnt, (double)good_cnt/all_cnt*100);
	printf("Novel transcripts with %s: %d (%.2lf%% of novel transcripts)\n",threshold_string, novel_good_t,(double)novel_good_t/novel_t*100);
	printf("Reference transcripts with %s: %d (%.2lf%% of reference transcripts)\n",threshold_string, good_ref_cnt,(double)good_ref_cnt/ref_cnt*100);

//check genes according to tracking file:
	max_gene_cnt = 0;
	for (i = 0; i < all_cnt; ++i) {
		if (evidence_table[i].gene[0] != '\0') {
			max_gene_cnt++;
		}
	}
	if ((tracking_gene_table = calloc(max_gene_cnt, sizeof(gene_evidence))) == NULL) {
                printf("Memory allocation error in main()\n");
                exit(0);
        }	
	for (i = 0, j=0; i < all_cnt; ++i) {
		if (evidence_table[i].gene[0] != '\0') {
			strcpy(tracking_gene_table[j].gene,evidence_table[i].gene);
			tracking_gene_table[j].found_in_rna = evidence_table[i].found_in_rna;
			j++;
		}
	}
	if (j != max_gene_cnt) {
		printf("Error in counting genes\n");
		exit(0);
	}
	qsort(tracking_gene_table, max_gene_cnt, sizeof(gene_evidence), sort_gene_alpha);
	gene_nodup_cnt = 1;
	for (i = 1; i < max_gene_cnt; ++i) {
		if (strcmp(tracking_gene_table[i].gene, tracking_gene_table[i-1].gene) != 0)
			gene_nodup_cnt++;
	}
	if (gene_nodup_cnt > 0) {
		if ((no_dup_tracking_gene_table = calloc(gene_nodup_cnt, sizeof(gene_evidence))) == NULL) {
        	        printf("Memory allocation error in main()\n");
               		exit(0);
        	}	
	}
	if (tracking_gene_table[0].gene[0] != '\0') {
		strcpy(no_dup_tracking_gene_table[0].gene,tracking_gene_table[0].gene);
		no_dup_tracking_gene_table[0].found_in_rna += tracking_gene_table[0].found_in_rna;
	}
	
	for (i = 1, j = 0; i < max_gene_cnt; ++i) {
		if (strcmp(tracking_gene_table[i].gene, tracking_gene_table[i-1].gene) == 0) {
			no_dup_tracking_gene_table[j].found_in_rna += tracking_gene_table[i].found_in_rna;
		} else {
			j++;
			no_dup_tracking_gene_table[j].found_in_rna += tracking_gene_table[i].found_in_rna;
			strcpy(no_dup_tracking_gene_table[j].gene,tracking_gene_table[i].gene);
		}	
	}
	for (i = 0; i < gene_nodup_cnt; ++i) {
		if (no_dup_tracking_gene_table[i].found_in_rna)
			good_gene_cnt++;
	}
	printf("Total genes in tracking file: %d\n\n",gene_nodup_cnt);

	//use good_transcripts to filter gtf, produce new gtf output
	if ((f = fopen(gtf_fn,"r")) == NULL) {
        	printf("Can't open file %s\n",gtf_fn);
                exit(0);
        }
	stem[0] = output_fn[0] = '\0';
       	strcpy(stem, gtf_fn);
	if ((pnt = strstr(stem, ".gtf")) == NULL) {
		printf("gtf fn needs to end .gtf\n");
		exit(0);
	}
	*pnt = '\0';
	sprintf(output_fn, "%s_%s.gtf", stem, (mode == FPKM) ? "FPKM" : "confhi");
	if ((g = fopen(output_fn,"w")) == NULL) {
                printf("Can't open file %s\n", output_fn);
                exit(0);
        }
	output_fn[0] = '\0';
	sprintf(output_fn, "%s_%s_table.csv", stem,(mode == FPKM) ? "FPKM" : "confhi");
	if ((table_output = fopen(output_fn,"w")) == NULL) {
                printf("Can't open file %s\n", output_fn);
                exit(0);
        }
//this prints new gtf and also saves info on whether transcript was in original (unfiltered) gtf	
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (is_good_transcript(line, &evidence_table, all_cnt)) {
			fprintf(g, "%s", line);
		}	
	}
	fprintf(table_output,"Transcript,Gene,over threshold,is reference transcript,present in unfiltered gtf\n");
	for (i=0; i < all_cnt; i++) {
		fprintf(table_output, "%s,%s,%c,%c,%c\n",evidence_table[i].transcript,evidence_table[i].gene,(evidence_table[i].found_in_rna) ? 'Y':'N', (evidence_table[i].is_ref_transcript) ? 'Y':'N', (evidence_table[i].in_original_gtf) ? 'Y':'N');
	}
	fclose(table_output);
//get more stats
	max_gene_cnt = 0, good_gene_cnt = 0, gene_nodup_cnt = 0;
	for (i = 0; i < all_cnt; ++i) {
		if(evidence_table[i].in_original_gtf) {
			unfilt_gtf_t++;
			if (evidence_table[i].is_ref_transcript) {
				final_ref_cnt++;
				if (evidence_table[i].found_in_rna) {
					final_good_ref_cnt++;
				}
			}
			if (evidence_table[i].found_in_rna) {
				final_t++;
			}
			if (evidence_table[i].gene[0] != '\0') {
				max_gene_cnt++;
			}
		}
	}
	printf("Unfiltered (input) gtf information:\n");
	printf("Total transcripts: %d\n",unfilt_gtf_t);
	printf("Novel transcripts: %d (%.2lf%% of total transcripts)\n", (unfilt_gtf_t - final_ref_cnt),(double)(unfilt_gtf_t - final_ref_cnt)/unfilt_gtf_t*100);

	if ((gene_table = calloc(max_gene_cnt, sizeof(gene_evidence))) == NULL) {
                printf("Memory allocation error in main()\n");
                exit(0);
        }	
	for (i = 0, j=0; i < all_cnt; ++i) {
		if (evidence_table[i].in_original_gtf && evidence_table[i].gene[0] != '\0') {
			strcpy(gene_table[j].gene,evidence_table[i].gene);
			gene_table[j].found_in_rna = evidence_table[i].found_in_rna;
			j++;
		}
	}
	if (j != max_gene_cnt) {
		printf("Error in processing gene information\n");
		exit(0);
	}
	qsort(gene_table, max_gene_cnt, sizeof(gene_evidence), sort_gene_alpha);
	
	gene_nodup_cnt = 1;
	for (i = 1; i < max_gene_cnt; ++i) {
		if (strcmp(gene_table[i].gene, gene_table[i-1].gene) != 0)
			gene_nodup_cnt++;
	}
	printf("Genes: %d\n\n", gene_nodup_cnt);
	if (gene_nodup_cnt > 0) {
		if ((no_dup_gene_table = calloc(gene_nodup_cnt, sizeof(gene_evidence))) == NULL) {
        	        printf("Memory allocation error in main()\n");
               		exit(0);
        	}	
	}
	if (gene_table[0].gene[0] != '\0') {
		strcpy(no_dup_gene_table[0].gene,gene_table[0].gene);
		no_dup_gene_table[0].found_in_rna += gene_table[0].found_in_rna;
	}
	
	for (i = 1, j = 0; i < max_gene_cnt; ++i) {
		if (strcmp(gene_table[i].gene, gene_table[i-1].gene) == 0) {
			no_dup_gene_table[j].found_in_rna += gene_table[i].found_in_rna;
		} else {
			j++;
			no_dup_gene_table[j].found_in_rna += gene_table[i].found_in_rna;
			strcpy(no_dup_gene_table[j].gene,gene_table[i].gene);
		}	
	}
	
	for (i = 0; i < gene_nodup_cnt; ++i) {
		if (no_dup_gene_table[i].found_in_rna)
			good_gene_cnt++;
	}
	printf("Filtered (output) gtf information:\n");
	printf("Total transcripts: %d\n", final_t);
	printf("Novel transcripts: %d (%.2lf%% of total transcripts)\n", (final_t-final_good_ref_cnt),(double)(final_t-final_good_ref_cnt)/final_t*100);
	printf("Percentage coverage of reference transcripts: %.2lf\n",(double)final_good_ref_cnt/final_ref_cnt*100);
	
	printf("Genes: %d\n",good_gene_cnt);
	printf("Percentage coverage of reference genes: %.2lf\n",(double)good_gene_cnt/gene_nodup_cnt*100);
//clean-up
	if (evidence_table != NULL)
		free(evidence_table);
	if (gene_table != NULL)
		free(gene_table);
	if (no_dup_gene_table != NULL)
		free(no_dup_gene_table);
	if (tracking_gene_table != NULL)
		free(tracking_gene_table);
	if (no_dup_tracking_gene_table != NULL)
		free(no_dup_tracking_gene_table);
	
	fclose(f);
	fclose(g);
	exit(0);
}
/************************************************************************/

/************************************************************************/
int is_good_transcript(char *line, evidence **transcripts, int cnt)
{
	char *pnt = NULL, *ptr = NULL, copy[MAXSTR], transcript[100];
	int i, transcript_num = 0;

	copy[0] = transcript[0] = '\0';
	strcpy(copy, line);
	//find transcript_id
	if ((ptr = strstr(copy, "transcript_id")) == NULL) {
		printf("This line has no transcript id - Abort\n%s",line);
		exit(0);
	}
	ptr += 15;
	pnt = move_ptr(ptr,'"');
	strcpy(transcript,ptr);

	for (i = 0; i < strlen(transcript); ++i) {
		if (!isdigit(*ptr)) {
			ptr += 1;
		} else {
			transcript_num = atoi(ptr);
			break;
		}
	}
	*pnt = '"';
	//check if transcript id is in list
	if (strcmp((*transcripts)[transcript_num-1].transcript,transcript) == 0) {
		//record that this transcript was in the original gtf
		(*transcripts)[transcript_num-1].in_original_gtf = 1;
		//print to new gtf if there is evidence
		if ((*transcripts)[transcript_num-1].found_in_rna) {
			return 1;
		}
	} else { //if transcript wasn't easily found, search for it by stepping through array
		for (i = 0; i < cnt; i++) {
			if (strcmp((*transcripts)[i].transcript,transcript) == 0) {
				(*transcripts)[i].in_original_gtf = 1;
				if ((*transcripts)[transcript_num-1].found_in_rna) {
					return 1;
				}
			}
		}
	}
	return 0;
}
/************************************************************************/

/************************************************************************/
void print_help() 
{
	printf("Required Input:\n\t-g\tcuffcompare output .combined.gtf file\n");
	printf("\t-t\tCuffcompare output tracking file\n");
	printf("Optional Input:\n\t-f\tthreshold value for output transcript selection\n");
	printf("\t\t\t(default = 0.0)\n");
	printf("\t-c\tset to select transcripts based on conf_hi value\n");
	printf("\t\t\t(default = select transcripts based on FPKM value)\n");
	printf("\t-h\tprint help\n");
	return;
}
/************************************************************************/

/************************************************************************/
evidence *get_good_transcripts(FILE *f, double threshold, int total, mode_type mode)
{	
	char line[MAXSTR], *pnt = NULL, *ptr = NULL, transcript[100];
	char copy[MAXSTR], class_code, gene[50];
	double fpkm = 0.0, conf_hi = 0.0;
	int good_one, line_num = 0;
	evidence *evidence_table = NULL;

	transcript[0] = gene[0] = '\0';
	if ((evidence_table = calloc(total, sizeof(evidence))) == NULL) {
                printf("Memory allocation error in get_good_transcripts()\n");
                exit(0);
        }
	//printf("Transcripts with %s > %lf:\n", (mode == FPKM) ? "FPKM":"conf_hi", threshold);
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		good_one = 0;
		gene[0] = '\0';
		strcpy(copy,line);
		ptr = copy;
		//get transcript
		pnt = move_ptr(ptr, '\t');
		strcpy(transcript, ptr);
		strcpy(evidence_table[line_num].transcript,ptr);
		*pnt = '\t';
		//get ref_gene
		ptr = skip_fields(ptr, '\t', 2);
		if (ptr[0] != '-') { 
			pnt = move_ptr(ptr, '|');
			strcpy(gene, ptr);
			*pnt = '|';
			strcpy(evidence_table[line_num].gene,gene);
		}
		//get class_code
		ptr = skip_fields(ptr, '\t', 1);
		pnt = move_ptr(ptr, '\t');
		class_code = ptr[0];
		if (class_code == '=' || class_code == 'c') {
			evidence_table[line_num].is_ref_transcript = 1;
		} 
		*pnt = '\t';
		ptr = pnt + 1;

		//get FPKM
		ptr = skip_fields(ptr, '|', 3);
		pnt = move_ptr(ptr, '|');
		fpkm = atof(ptr);
		*pnt = '|';
		ptr = pnt + 1;

		//get conf_hi
		ptr = skip_fields(ptr, '|', 1);
		pnt = move_ptr(ptr, '|');
		conf_hi = atof(ptr);
		*pnt = '|';

		//save good ones
		if (mode == FPKM) {
			if (isgreater(fpkm, threshold)) {
				good_one = 1;
			}
		} else if (mode == CONF_HI) {	
			if (isgreater(conf_hi, threshold)) { 
				good_one = 1;
			}
		}

		if (good_one) {
			//print log
		//	printf("%s",line);
		//	fflush(stdout);	
			//save transcript
			evidence_table[line_num].found_in_rna = 1;
		} //else if (!(class_code == '=' || class_code == 'c') && mode == CONF_HI) {
		//	printf("This one is novel but no good? class_code %c\n%s", class_code, line);  
		//}
		line_num++;
	}
	printf("\n");
	return(evidence_table);
}
/************************************************************************/

/************************************************************************/
int sort_gene_alpha(const void *a, const void *b) 
{
	gene_evidence *e1 = (gene_evidence *)a, *e2 = (gene_evidence *)b;
    	return strcmp(e1->gene, e2->gene);
}
/************************************************************************/

/************************************************************************/
int count_lines(FILE *f)
{
        int cnt = 0;
        char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
        	cnt++;
        }
        rewind(f);

        return(cnt);
}
/************************************************************************/

/************************************************************************/
char *move_ptr(char *str,char delim)
{
	char *tmp = NULL;

	if((tmp = strchr(str,delim)) != NULL) *tmp = '\0';

	return(tmp);
}
/************************************************************************/

/************************************************************************/
char *skip_fields(char *str, char delim, int num)
{
        char *tmp = NULL, *ptr;
        int i;

        ptr = str;
        for (i = 0; i < num; ++i) {
                if ((tmp = strchr(ptr,delim)) != NULL && *(tmp + 1) != '\0') {
                        ptr = tmp + 1;
                } else {
                        printf("Can't skip %d fields\n", num);
                }
        }

        return(ptr);
}
/************************************************************************/
