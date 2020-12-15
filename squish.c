#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>
#define QFRAG 3
#define MAXSTR 51200
#define FASTALEN 70
#define MINLEN 8
#define MAXTHREADS 8
#define CATLEN 300
#define SEPARATOR "WWWWW"

typedef struct param_info param_info;
typedef struct find_info find_info;
typedef struct entry entry;
typedef struct seq_node seq_node;
typedef struct qfrag_node qfrag_node;

struct find_info {
	int seq_num;
	int pos;
};
struct entry {
        char *header;
        char *seq;
        int len;
        int collapse;
	int reps;
	find_info find;
	char collapse_info[100];
	long hash_start;
};
struct qfrag_node {
	long *loc;
	long cnt;
};
struct param_info {
	int num;
	int max;
	entry *no_dup;
	int no_dup_cnt;
	int *skip_array;
	int skip_array_cnt;
	qfrag_node *hash_table;
};

char aa[] = "ACDEFGHIKLMNPQRSTVWXY";

int count_seq(FILE *);
entry *get_seq(FILE *, int);
int hash_this(char *);
qfrag_node *create_hash_table(entry **, int, int, long);
int sort_alpha(const void *, const void *);
int sort_by_len(const void *, const void *);
int *get_qfrags(char *, int *, int *, int *); 
find_info find_last(qfrag_node *, int, entry *, int, int *, int, int, int);
long *and_loc(long *, qfrag_node, qfrag_node, int, int);
long find_loc(long, qfrag_node, long);
qfrag_node *copy_qfrag_node(qfrag_node *, long);
void remove_duplicates(entry **, int *);
void merge_headers(char *, char **);
void add_meta(char *h1, char **h2, int st, int en);
void print_fasta(FILE *, entry);
int *make_skip_array(entry *, int, int *);
void substring_search(void *);
void free_entries(entry **, int); 
void merge_entries(entry **, int *, entry *, int);
void make_headers_unique(entry **, int, int);
void print_usage();
/*
squish.c: merge protein fasta files, remove duplicate and use a 'hash table' index to search for and compact redundancy
	specify protein fasta files with -d; number of threads with -t; output file name with -o. Default output
	file name is 'output.fasta'.
*/
/**********************************************************************************************************/
void main(int argc,char **argv)
{
	FILE *f, *g;
	entry *entries = NULL, *next_entries = NULL, artificial;
	find_info find;
	qfrag_node *hash_table;
	int hash = 0, hash_num = (int)(pow(strlen(aa),QFRAG)), cnt = 0, next_cnt = 0, final_cnt = 0, ncqf_cnt = 0, i, j, k, index, *ncqf = NULL, last = -1, mod = 0, no_dup_cnt = 0, removed_cnt = 0, skip = 0, *skip_array = NULL, skip_array_cnt = 0, thread_num = 1, input_cnt = 0, no_dup2_cnt = 0, cat_cnt = 0;
	long hash_len = 0, test_len = 0;	
	char **input, **realloc_tmp, output[100], *pnt = NULL, c, tmp_header[MAXSTR], artificial_seq[MAXSTR], artificial_header[MAXSTR];
	pthread_t *threads = NULL;
	param_info *params = NULL;

	output[0] = '\0';
	strcpy(output, "output.fasta");
	while ((c = getopt (argc, argv, "t:d:o:h")) != -1) {
                switch (c) {
                        case 't':
                                thread_num = atoi(optarg);
                                printf("Number of threads specified: %d\n", thread_num);
				if (thread_num > MAXTHREADS) {
                			printf("%d is too many threads. MAXTHREADS is currently %d\n",thread_num, MAXTHREADS);
                			exit(0);
				}
                                break;
			case 'd':
				if (!input_cnt) {
        				if ((input = calloc(1,sizeof(char *))) == NULL) {
        	        			printf("Memory allocation error in main()\n");
        	        			exit(0);
        				}
				}
				else {
                                	if ((realloc_tmp = realloc(input, (input_cnt + 1) * sizeof(char *))) == NULL) {
                                        	printf("memory allocation error in main()\n");
                                                exit(0);
                                        }
                                        input = realloc_tmp;
                                }
        			if ((input[input_cnt] = calloc(strlen(optarg) + 1,sizeof(char))) == NULL) {
        			        printf("Memory allocation error in main()\n");
        			        exit(0);
        			}
				strcpy(input[input_cnt], optarg);
				input_cnt ++;
				break;
			case 'o':
				output[0] = '\0';
				strcpy(output, optarg);
				break;
			case 'h':
				print_usage();
				exit(0);
                        default: /*'?'*/
                                if (optopt == 't')
                                        printf ("With option -t, specify the number of threads\n\n");
				else if (optopt == 'd')
                                        printf ("With option -d, specify a protein fasta database\n\n");
				else if (optopt == 'o')
                                        printf ("With option -o, specify an output file name\n\n");
                                else
                                        printf("Unknown option '-%c'.\n\n", optopt);
				print_usage();
                                exit(0);
                }
        }
	if (optind != argc) {
                printf("Check your command. All specifications following the program name should be via options (-d, -o, -t)\n\n");
		print_usage();
                exit(0);
        }

	if (!input_cnt) {
		printf("Please specify one or more protein fasta databases using -d option\n\n");
		print_usage();
		exit(0);
	}
	for (i = 0; i < input_cnt; ++i) { 
		printf("input file %d: %s\n", i + 1, input[i]);
		if (strstr(input[i], ".fa") == NULL) {
			printf("input databases must be .fasta or .fa format\n\n");
			print_usage();
			exit(0);
		}
	}
	if ((g = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n", output);
                exit(0);
        }
	if ((f = fopen(input[0],"r")) == NULL) {
               	printf("Can't open file %s\n",input[0]);
               	exit(0);
        }
	
	cnt = count_seq(f);
        entries = get_seq(f, cnt);
	fclose(f);
        printf("%d entries in %s\n", cnt, input[0]);
	
	make_headers_unique(&entries, cnt, 1);
	qsort(entries, cnt, sizeof(entry), sort_alpha);
        qsort(entries, cnt, sizeof(entry), sort_by_len);
//delete duplicates
        remove_duplicates(&entries, &cnt);
	
	printf("%d entries left in %s after remove_duplicates()\n", cnt, input[0]);
//get data from next file, remove duplicates, then merge with previous data and remove duplicates
	for (i = 1; i < input_cnt; ++i) {
		if ((f = fopen(input[i],"r")) == NULL) {
                	printf("Can't open file %s\n",input[i]);
                	exit(0);
        	}
        	next_cnt = count_seq(f);
        	next_entries = get_seq(f, next_cnt);
        	fclose(f);

        	printf("%d entries in %s\n", next_cnt, input[i]);
		make_headers_unique(&next_entries, next_cnt, i+1);
		qsort(next_entries, next_cnt, sizeof(entry), sort_alpha);
        	qsort(next_entries, next_cnt, sizeof(entry), sort_by_len);
      
	  	remove_duplicates(&next_entries, &next_cnt);
		printf("%d entries left in %s after remove_duplicates()\n", next_cnt, input[i]);
	
		merge_entries(&entries, &cnt, next_entries, next_cnt);
        	printf("After merge, %d entries stored from %d files\n", cnt, i+1);
		
		free_entries(&next_entries, next_cnt);
		next_entries = NULL;
	}
//add removed duplicate metadata to headers
	for (i = 0; i < cnt; ++i) {
		if (entries[i].reps) {
			tmp_header[0] = '\0';
			sprintf(tmp_header, "%s had %d exact replicates", entries[i].header, entries[i].reps);
			free(entries[i].header);
        		if ((entries[i].header = calloc(strlen(tmp_header) + 1,sizeof(char))) == NULL) {
                		printf("Memory allocation error in remove_duplicates()\n");
                		exit(0);
        		}
			strcpy(entries[i].header, tmp_header);
		}	
	}
//squish
	for (i = 0; i < cnt; ++i) {
		hash_len += entries[i].len - (QFRAG - 1);
		if (hash_len > (LONG_MAX - 100000)) {
                       	printf("cumulative sequence length getting close to LONG_MAX - change everything to long long!\n");
                       	exit(0);
		}
	}
	printf("hash_len = %ld\n", hash_len);
	
	hash_table = create_hash_table(&entries, cnt, hash_num, hash_len); 

//make skip_array
	skip_array = make_skip_array(entries, cnt, &skip_array_cnt);

//chop search seq into non-contiguous Qfrags - make array. - if strlen(seq)%QFRAG!=0, add final qfrag left-shifted
	if (skip_array_cnt) {
        	if ((threads = calloc(thread_num,sizeof(pthread_t))) == NULL) {
        	        printf("Memory allocation error in main()\n");
        	        exit(0);
        	}
        	if ((params = calloc(thread_num,sizeof(param_info))) == NULL) {
        	        printf("Memory allocation error in main()\n");
        	        exit(0);
        	}
		
		for (i = 0; i < thread_num; ++i) {
			params[i].num = i;
			params[i].max = thread_num;
			params[i].no_dup = &entries[0];
			params[i].no_dup_cnt = cnt;
			params[i].skip_array = &skip_array[0];
			params[i].skip_array_cnt = skip_array_cnt;
			params[i].hash_table = &hash_table[0];
			pthread_create(&threads[i], NULL, (void *)substring_search, (void *)&params[i]);
		}
		for (i = 0; i < thread_num; ++i)
			pthread_join(threads[i], NULL);	
	}
        
	for (i = 0; i < cnt; ++i) {
		if (entries[i].collapse) {
			add_meta(entries[i].header, &entries[entries[i].find.seq_num].header, entries[i].find.pos, (entries[i].find.pos + entries[i].len - 1));
		}
	}
	artificial.header = NULL;
	artificial.seq = NULL;
	artificial.len = 0;
	artificial_seq[0] = '\0';
	artificial_header[0] = '\0';
	cat_cnt = 0;
	for (i = 0; i < cnt; ++i) {
		if (!entries[i].collapse) {
		//concatenate short sequences
			if (entries[i].len < CATLEN + 1) {
				cat_cnt++;
				if (cat_cnt == 1) {
					strcpy(artificial_seq, entries[i].seq);
					strcpy(artificial_header, entries[i].header); 
				} else {
					strcat(artificial_seq, SEPARATOR);
					strcat(artificial_seq, entries[i].seq);
				}
				if (strlen(artificial_seq) >= (CATLEN * 2) || i == cnt - 1 || entries[i+1].len > CATLEN) {
        				if ((artificial.seq = calloc(strlen(artificial_seq) + 1,sizeof(char))) == NULL) {
        	        			printf("Memory allocation error in main()\n");
        	        			exit(0);
        				}
        				if ((artificial.header = calloc(strlen(artificial_header) + 25,sizeof(char))) == NULL) {
        	        			printf("Memory allocation error in main()\n");
        	        			exit(0);
        				}
					strcpy(artificial.seq, artificial_seq);
					sprintf(artificial.header, "%s + %d concatenations", artificial_header, cat_cnt);
					
					print_fasta(g, artificial);
					final_cnt++;

					//reset artificial entry and temp header and seq
					if (artificial.header != NULL)
						free(artificial.header);
					if (artificial.seq != NULL)
						free(artificial.seq);
					artificial.header = NULL;
					artificial.seq = NULL;
					artificial_header[0] = '\0';
					artificial_seq[0] = '\0';
					cat_cnt = 0;
				}
			} else {	
				print_fasta(g, entries[i]);
				final_cnt++;
			}
		}
		else {
			removed_cnt++;
		}
	}
	printf("Removed %d sequences because they were wholly contained within another sequence.\nShort sequences were then concatenating for printing.\nThus, %s contains %d protein sequences\n", removed_cnt, output, final_cnt);
	free_entries(&entries, cnt);
	for (i = 0; i < input_cnt; ++i)
		if (input[i] != NULL)
			free(input[i]);
	if (input != NULL)
		free(input);
	for (i = 0; i < hash_num; ++i) { //for each qfrag_node
		if (hash_table[i].cnt>0) {
			if (hash_table[i].loc != NULL)
				free(hash_table[i].loc);
		}
	}
	if (hash_table != NULL)
		free(hash_table);
	if (threads != NULL)
		free(threads);
	if (params != NULL)
		free(params);
	if (skip_array != NULL)
		free(skip_array);
	fclose(g);
	exit(0);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void print_usage()
{
	printf("Required input:\n\t-d\tspecify one protein fasta database\n\t\t(this option should be repeated to add more databases;\n\t\te.g. –d database1.fasta –d database2.fasta)\n");
	printf("Optional:\n\t-t\tSpecify number of threads to use [default = 1]\n");
	printf("\t-o\tSpecify name of output database file [default = \"output.fasta\"]\n");
	printf("\t-h\tprint usage help\n\n");
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void make_headers_unique(entry **entries, int cnt, int file_num) 
{
	int i, j;
	char *tmp, *ptr_from = NULL, *ptr_to = NULL, append[10];

	append[0] = '\0';
	sprintf(append, "_%d", file_num);
	for (i = 0; i < cnt; ++i) {
		ptr_from = (*entries)[i].header;
        	if ((tmp = calloc(strlen((*entries)[i].header) + 3,sizeof(char))) == NULL) {
        	        printf("Memory allocation error in make_headers_unique()\n");
        	        exit(0);
        	}
		ptr_to = tmp;
		while (*ptr_from != '\0' && *ptr_from != ' ') {
			*ptr_to++ = *ptr_from++;
		}
		strcat(tmp, append);
		ptr_to += strlen(append);
		while (*ptr_from != '\0') {
			*ptr_to++ = *ptr_from++;
		} 
		*ptr_to = '\0';
		if ((*entries)[i].header != NULL)
			free((*entries)[i].header);
		(*entries)[i].header = tmp;
	}

	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void substring_search(void *p)
{
	param_info *params = (param_info *)p;
	int ncqf_cnt = 0, i, j, k, *ncqf = NULL, last = -1, mod = 0, skip = 0, ctr = 0;


	printf("Thread %d up and running\n", params->num);
	skip = params->skip_array[0];
        for (i = 0, j=0; i < params->no_dup_cnt - 1; ++i) {
        	if (ctr == params->num) {
			if (i%100000 == 0) {
                                printf("Thread %d searching for seq %d\n", params->num, i+1);
                                fflush(stdout);
                        }
                        ncqf = get_qfrags(params->no_dup[i].seq, &ncqf_cnt, &last, &mod);
			
                        if (params->hash_table[ncqf[0]].cnt) {
				while (i > skip) {
					j++;
                                        if (j < params->skip_array_cnt) {
                                                skip = params->skip_array[j];
                                        }
                                        else {
                                                if (ncqf != NULL)
                                                        free(ncqf);
                                                break;
                                        }
                                }
                                if (i == skip) {
                                        j++;
                                        if (j < params->skip_array_cnt) {
                                                skip = params->skip_array[j];
                                        }
                                        else {
                                                if (ncqf != NULL)
                                                        free(ncqf);
						ncqf = NULL;
                                                break;
                                        }
                                }
				
                                params->no_dup[i].find = find_last(params->hash_table, skip, params->no_dup, params->no_dup_cnt, ncqf, ncqf_cnt, last, mod);
                                if (params->no_dup[i].find.seq_num != -1) {
                                        params->no_dup[i].collapse = 1;
                                }
//reset
                        }
                        ncqf_cnt = 0;
                        if (ncqf != NULL)
                                free(ncqf);
			ncqf = NULL;
                        last = -1;
                        mod = 0;
                }
		ctr++;
		if (ctr == params->max)
			ctr = 0;
	}
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int *make_skip_array(entry *entries, int cnt, int *skip_array_cnt)
{
	int *tmp = NULL, i, len = entries[0].len, j=0;
	*skip_array_cnt = 0;

	for (i = 1; i < cnt; i++) {
		if (entries[i].len > len) {
			*skip_array_cnt+=1;
			len = entries[i].len;
		}
	}
	if (*skip_array_cnt) {
        	if ((tmp = calloc(*skip_array_cnt,sizeof(int))) == NULL) {
        	        printf("Memory allocation error in make_skip_array()\n");
        	        exit(0);
        	}
		len = entries[0].len;
		for (i = 0; i < cnt; i++) {
			if (entries[i].len > len) {
				tmp[j++] = i;
				len = entries[i].len;
			}	
		}
	}
	return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void print_fasta(FILE *f, entry e)
{
        int len = (int)strlen(e.seq), i = 0;
        char copy[FASTALEN + 2];

        fprintf(f, "%s\n", e.header);
        while (i < len) {
                copy[0] = '\0';
                strncpy(copy, &(e.seq[i]), FASTALEN);
                copy[FASTALEN] = '\0';
                fprintf(f, "%s\n", copy);
                i += FASTALEN;
        }
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void add_meta(char *h1, char **h2, int st, int en)
{
        char tmp[MAXSTR];

        tmp[0] = '\0';
        if ((strlen(h1) + strlen(*h2) +50) > MAXSTR) {
                printf("Increase MAXSTR to accommodate new header\n");
                exit(0);
        }

        sprintf(tmp, "%s (%s@%d->%d);", *h2, h1, st, en);
        if (*h2 != NULL)
                free(*h2);
        *h2 = NULL;
        if ((*h2 = (char *)calloc((strlen(tmp) +1),sizeof(char))) == NULL) {
                printf("Memory allocation error in merge_headers()\n");
                exit(0);
        }
        strcpy(*h2, tmp);

        return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int hash_this(char *seq)
{
	extern char aa[]; 
	char *ptr = NULL;
	int i,j, tmp = 0, index, len = strlen(aa), seq_len = strlen(seq);

	if (seq_len > QFRAG) {
		printf("code is too long: %s\n", seq);
		exit(0);	
	}	
	for (i = seq_len-1, j=0; i>-1; --i, j++) {
		if ((ptr = strchr(aa, seq[i])) != NULL) {
			index = len - strlen(ptr);
			tmp +=	(index * (int)pow(len, j));
		}
		else {
			printf("disallowed amino acid %c in seq %s\n", seq[i], seq);
			exit(0);
		}
	}
	
	return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int count_seq(FILE *f)
{
        int cnt = 0;
        char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>')
                        cnt++;
        }
        rewind(f);

        return(cnt);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
entry *get_seq(FILE *f, int cnt)
{
        char line[MAXSTR], seq[MAXSTR];
        entry *tmp = NULL;
        int first = 1, i, len = 0;

        if ((tmp = (entry *)calloc(cnt,sizeof(entry))) == NULL) {
                printf("Memory allocation error in get_seq()\n");
                exit(0);
        }
	i=0;
        seq[0] = '\0';
        while(fgets(line,MAXSTR - 1,f) != NULL) {
                while (line[(int)strlen(line) - 1] == '\n' || line[(int)strlen(line) - 1] == '\r')
                	line[(int)strlen(line) - 1] = '\0';
                if (line[0] == '>') {
                        if (first) {
                                if ((tmp[i].header = calloc(((int)strlen(line) + 1),sizeof(char))) == NULL) {
                                        printf("Memory allocation error in get_seq()\n");
                                        exit(0);
                                }
                                strcpy(tmp[i].header, line);
                                first = 0;
                        }
                        else {
				if ((tmp[i].seq = calloc((strlen(seq) + 1),sizeof(char))) == NULL) {
                                        printf("Memory allocation error in get_seq()\n");
                                        exit(0);
                                }
                                strcpy(tmp[i].seq, seq);
                                tmp[i].len = strlen(seq);
                                seq[0] = '\0';
                                len = 0;
                                i++;
                                if ((tmp[i].header = calloc((strlen(line) + 1),sizeof(char))) == NULL) {
                                        printf("Memory allocation error in get_seq()\n");
                                        exit(0);
                                }
                                strcpy(tmp[i].header, line);
                        }
                }
                else {
                        len += strlen(line);
                        if (len < MAXSTR) {
                                strcat(seq, line);
			}
			else {
                                printf("sequence len %d exceeds max length of %d - increase MAXSTR\n", len, MAXSTR);
                                exit(0);
                        }
                }
        }
        if ((tmp[i].seq = calloc((int)(strlen(seq) + 1),sizeof(char))) == NULL) {
                printf("Memory allocation error in get_seq()\n");
                exit(0);
        }
        strcpy(tmp[i].seq, seq);
        tmp[i].len = strlen(seq);
        i++;
        if (i != cnt)
                printf("What the?\n");

        return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int sort_by_len(const void *a, const void *b)
{
        entry *e1 = (entry *)a, *e2 = (entry *)b;

	return (e1->len - e2->len);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int sort_alpha(const void *a, const void *b) {
        entry *e1 = (entry *)a, *e2 = (entry *)b;
    return strcmp(e1->seq, e2->seq);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
qfrag_node *create_hash_table(entry **entries, int entries_cnt, int node_cnt, long hash_len)
{
	char qf[QFRAG +1];
	qfrag_node *hash_table = NULL;
	int pos, seq_num, i, hash;
	long loc = 0, *tmp = NULL;	

	printf("creating hash table...\n");
        if ((hash_table = (qfrag_node *)calloc(node_cnt,sizeof(qfrag_node))) == NULL) {
                printf("Memory allocation error in create_hash_table()\n");
                exit(0);
        }
	for (i=0; i<node_cnt; ++i) {
		if ((hash_table[i].loc = (long *)calloc(hash_len,sizeof(long))) == NULL) {
                	printf("Memory allocation error in create_hash_table()\n");
                	exit(0);
        	}
	}
	for (seq_num = 0; seq_num < entries_cnt; seq_num++) {
		if (seq_num%100000 == 0) {
			printf("hashing sequence %d of %d\n", seq_num+1, entries_cnt);
			fflush(stdout);
		}
		(*entries)[seq_num].hash_start = loc;
		for (pos = 0; pos < ((*entries)[seq_num].len - QFRAG +1); pos++) {
			strncpy(qf, &(*entries)[seq_num].seq[pos], QFRAG);
			qf[QFRAG] = '\0';
			hash = hash_this(qf);
			hash_table[hash].loc[hash_table[hash].cnt++] = loc++;		
		}
	}
	if (loc !=hash_len)
		printf("eh?? loc = %ld, hash_len = %ld\n", loc, hash_len);
//gather wasted memory
	for (i=0; i<node_cnt; ++i) {
		if (hash_table[i].cnt) {
			if ((tmp = (long *)calloc(hash_table[i].cnt,sizeof(long))) == NULL) {
                		printf("Memory allocation error in create_hash_table()\n");
                		exit(0);
        		}
			memcpy(tmp, hash_table[i].loc, hash_table[i].cnt * sizeof(long));
			if (hash_table[i].loc != NULL)
				free(hash_table[i].loc);
			hash_table[i].loc = tmp;
		}
		else {
			if (hash_table[i].loc != NULL)
				free(hash_table[i].loc);
		}
	}
	
	return hash_table;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int *get_qfrags(char *seq, int *cnt, int *last, int *mod)
{ 
	int *tmp = NULL;
	int i;
	char code[QFRAG+1];

	code[0] = '\0';
	if ((tmp = (int *)calloc((int)(strlen(seq)/3),sizeof(int))) == NULL) {
        	printf("Memory allocation error in get_qfrags()\n");
                exit(0);
	}

	for (i = 0; i< ((int)strlen(seq) - QFRAG + 1); i+=QFRAG) {
		strncpy(code, &seq[i], QFRAG);
		code[QFRAG] = '\0';
		tmp[*cnt] = hash_this(code);
		code[0] = '\0';
		*cnt += 1;
	}
	
	if ((*mod = (int)strlen(seq) % QFRAG) != 0) {
		strncpy(code, &seq[strlen(seq) - QFRAG], QFRAG);
		*last = hash_this(code); 
	} 

	return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
find_info find_last(qfrag_node *hash_table, int skip, entry *entries, int entries_cnt, int *ncqf, int cnt, int last, int mod)
{ 
	find_info find;
	int j, k, iteration, seq_len = cnt * QFRAG + mod/*, len = entries[skip].len*/;
	long *tmp_loc = NULL, loc_match = 0, i, tmp = -1;
	qfrag_node *copy = NULL;

	find.seq_num = -1;
	find.pos = -1;

//get qfrag_node for hash. Trim successively by copying nodes only if they match relevant loc for next hash.
	if (hash_table[ncqf[0]].cnt) {
		copy = copy_qfrag_node(&hash_table[ncqf[0]], entries[skip].hash_start);
		for (iteration = 1; iteration < cnt; iteration++) {
			tmp_loc = and_loc(&loc_match, *copy, hash_table[ncqf[iteration]], iteration, 0); 
			if (copy->loc != NULL)
				free(copy->loc);
			copy->loc = tmp_loc;
			copy->cnt = loc_match;	
			if (!copy->cnt) 
				break;
		}
		if (copy->cnt && last != -1) {
			tmp_loc = and_loc(&loc_match, *copy, hash_table[last], iteration, mod); 
			if (copy->loc != NULL)
				free(copy->loc);
			copy->loc = tmp_loc;
			copy->cnt = loc_match;	
		}
//curate for locs that bridge entries and get find_info
		for (i = 0, j = skip; i < copy->cnt; i++) { //for each loc left in copy
                	while (j != entries_cnt && (entries[j].hash_start < (copy->loc[i] + 1)))
                        	j++; //skip entries till you get the one with the loc in it
                	j--;
                	if ((copy->loc[i] + seq_len) < (entries[j].hash_start + entries[j].len + 1)) {
				tmp = copy->loc[i];
                        	find.seq_num = j;
                	}	
        	}
        	if (tmp != -1)
                	find.pos = (int)(tmp - entries[find.seq_num].hash_start);
//free copy
		if (copy->loc != NULL) 
			free(copy->loc);
		if (copy != NULL)
			free(copy);
	}

	return find;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
long *and_loc(long *match, qfrag_node old, qfrag_node new, int iteration, int mod)
{
	long start = 0, find = -1, tmp[old.cnt], *found = NULL, sought = 0, i;

	*match = 0;
	
	for (i = 0; i < old.cnt; ++i) {
		if (mod > 0) {
			sought = old.loc[i] + (iteration * QFRAG) - (QFRAG - mod);
		}
		else {
			sought = old.loc[i] + (iteration * QFRAG); 
		}
		if ((find = find_loc(sought, new, start)) != -1) {
			tmp[(*match)] = old.loc[i];
			*match+=1;
			start = find + 1;
			if (start == new.cnt) {
				break;
			}
		} 
	} 
	if (*match > 0) {
		if ((found = (long *)calloc(*match,sizeof(long))) == NULL) {
                	printf("Memory allocation error in and_loc()\n");
                	exit(0);
	        }	
		memcpy(found, &tmp, (*match) * sizeof(long));
	}
	
	return found;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
long find_loc(long sought, qfrag_node node, long start)
{
	long find = -1, i, l = start, r = node.cnt-1, mid = 0;

	if (l == r) {
		if (node.loc[l] == sought)
			find = l;
	}
	else {
		while ((r - l) != 1) {
			mid = (long)(((r-l)/2) + l);
			if (node.loc[mid] == sought) {
				find = mid;
				break;
			}
			else {
				if (node.loc[mid] < sought) {
					l = mid;
				}
				else {
					r = mid;
				}
			}
		}
		if (node.loc[l] == sought)
			find = l;
		if (node.loc[r] == sought)
			find = r;	
	}
	return find;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
qfrag_node *copy_qfrag_node(qfrag_node *template, long start)
{
	qfrag_node *tmp = NULL;
	long i, j, tmp_loc[template->cnt], exclude = 0;

	if ((tmp = (qfrag_node *)calloc(1,sizeof(qfrag_node))) == NULL) {
                printf("Memory allocation error in copy_qfrag_node()\n");
                exit(0);
        }
	for (i = 0; i < template->cnt; ++i) {
		if (template->loc[i] < start) {
			++exclude;
		}
	}
		
	tmp->cnt = template->cnt - exclude;
	if (tmp->cnt) {
		if ((tmp->loc = (long *)calloc(tmp->cnt,sizeof(long))) == NULL) {
                	printf("Memory allocation error in copy_qfrag_node()\n");
                	exit(0);
        	}
		for (i = exclude, j = 0; j < tmp->cnt; i++, j++) 
			tmp->loc[j] = template->loc[i];
	}

	return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void merge_entries(entry **entries, int *cnt, entry *next_entries, int next_cnt)
{
	entry *tmp = NULL;
	int i, total = *cnt + next_cnt;
	int j;
	
//copy all to new array
        if ((tmp = calloc(total,sizeof(entry))) == NULL) {
               	printf("Memory allocation error in merge_entries()\n");
               	exit(0);
        }
      	memcpy(&tmp[0], *entries, (*cnt * sizeof(entry)));
      	memcpy(&tmp[*cnt], next_entries, (next_cnt * sizeof(entry)));
	for (i = 0; i < *cnt; ++i) {
        	if ((tmp[i].header = calloc(strlen((*entries)[i].header) + 1,sizeof(char))) == NULL) {
                	printf("Memory allocation error in merge_entries()\n");
                	exit(0);
        	}
		strcpy(tmp[i].header, (*entries)[i].header);
        	if ((tmp[i].seq = calloc(strlen((*entries)[i].seq) + 1,sizeof(char))) == NULL) {
                	printf("Memory allocation error in merge_entries()\n");
                	exit(0);
        	}
		strcpy(tmp[i].seq, (*entries)[i].seq);
	}
	j = *cnt;	
	for (i = 0; i < next_cnt; ++i, ++j) {
        	if ((tmp[j].header = calloc(strlen(next_entries[i].header) + 1,sizeof(char))) == NULL) {
                	printf("Memory allocation error in merge_entries()\n");
                	exit(0);
        	}
		strcpy(tmp[j].header, next_entries[i].header);
        	if ((tmp[j].seq = calloc(strlen(next_entries[i].seq) + 1,sizeof(char))) == NULL) {
                	printf("Memory allocation error in merge_entries()\n");
                	exit(0);
        	}
		strcpy(tmp[j].seq, next_entries[i].seq);
	}
//sort and remove duplicates	
	qsort(tmp, total, sizeof(entry), sort_alpha);
        qsort(tmp, total, sizeof(entry), sort_by_len);
      	remove_duplicates(&tmp, &total);

	free_entries(entries, *cnt);
        *entries = tmp;
	*cnt = total;

	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void remove_duplicates(entry **entries, int *cnt)
{
        entry *tmp = NULL;
        int i, j = 0, k = 0, rep = 0, total = *cnt, index = 0;
	char current[MAXSTR], tmp_header[MAXSTR];

        for (i=0; i<total; ++i) {
		if (i == 0) {
			strcpy(current,(*entries)[0].seq);
			index = 0;
		}
		else {
                	if ((strcmp(current, (*entries)[i].seq)) == 0) {
                        	(*entries)[i].collapse = 1;
				if ((*entries)[i].reps)
					(*entries)[index].reps += (*entries)[i].reps;
				(*entries)[index].reps ++;
                        	j++;
                	}
			else {
				index = i;
				strcpy(current, (*entries)[i].seq);
			}
		}
        }
        *cnt = total - j;
		
//copy entries we want to keep
        if ((tmp = calloc((*cnt),sizeof(entry))) == NULL) {
                printf("Memory allocation error in remove_duplicates()\n");
                exit(0);
        }
        j=0;
        for (i=0; i<total; ++i) {
                if (!(*entries)[i].collapse) {
                        memcpy(&tmp[j], &(*entries)[i], sizeof(entry));
        		if ((tmp[j].header = calloc(strlen((*entries)[i].header) + 1,sizeof(char))) == NULL) {
                		printf("Memory allocation error in remove_duplicates()\n");
                		exit(0);
        		}
			strcpy(tmp[j].header, (*entries)[i].header);
        		if ((tmp[j].seq = calloc(strlen((*entries)[i].seq) + 1,sizeof(char))) == NULL) {
                		printf("Memory allocation error in remove_duplicates()\n");
                		exit(0);
        		}
			strcpy(tmp[j].seq, (*entries)[i].seq);
			j++;
                }
                else {
                        k++;
                }
        }
//free old entries:
	free_entries(entries, total);
        *entries = tmp;
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void free_entries(entry **entries, int cnt) 
{
	int i;

	for (i = 0; i < cnt; ++i) {	
        	if ((*entries)[i].seq != NULL)
                        free((*entries)[i].seq);
                if ((*entries)[i].header != NULL)
                	free((*entries)[i].header);
        }
        if ((*entries) != NULL)
                free(*entries);
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void merge_headers(char *h1, char **h2)
{
        int len = strlen(h1) + strlen(*h2) +2;
        char tmp[len];

        strcpy(tmp, h1);
        strcat(tmp, *h2);
        tmp[strlen(h1)] = '+';
        if (*h2 != NULL)
                free(*h2);
        *h2 = NULL;
        if ((*h2 = (char *)calloc(len,sizeof(char))) == NULL) {
                printf("Memory allocation error in merge_headers()\n");
                exit(0);
        }
        strcpy(*h2, tmp);

        return;
}
/**********************************************************************************************************/
