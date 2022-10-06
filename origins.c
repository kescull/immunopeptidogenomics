#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <curl/curl.h>
#include </home/kate/Downloads/cJSON/cJSON-master/cJSON.h>
// Download cJSON, update the filepath below to point to the correct folder and uncomment the line
/*#include </path/to/cJSON.h>*/

#define MAXSTR 51200
#define MAXSEQLEN 350000
#define MAXPROT 100000

typedef enum {
        INDEL,
        SUB,
        NO_CHANGE,
} Mut;

typedef enum {
        MINUS,
        PLUS,
        UNKNOWN,
} Strand;
typedef enum {
	NOT_FOUND,
	INTRON,
	UPSTREAM,
	DOWNSTREAM,
	UTR5,
	UTR3,
	TRANSLATED,
	NONCODING,
} Region;
 typedef enum{
       	EMBL_ON,
       	SIMPLE,
} Mode_type;

typedef struct tracking_info tracking_info;
typedef struct genome_info genome_info;
typedef struct seq_info seq_info;
typedef struct exon exon;
typedef struct vcf_mut_info vcf_mut_info;
typedef struct vcf_info vcf_info;
typedef struct contig_info contig_info;
typedef struct transcriptome_mutations transcriptome_mutations;
typedef struct mut_info mut_info;
typedef struct curl_data curl_data;
typedef struct genome_chunk genome_chunk;
typedef struct embl_info embl_info;
typedef struct transcriptome_info transcriptome_info;

struct tracking_info {
	char cuff_transcript[512];
        char nearest_ref[512];
        char class_code;
};
struct genome_chunk {
	int st;
	int en;
	int embl_st;
	int embl_en;
	int bases_from_start;
};
struct embl_info {
	char id[100];
	int id_num;
	int version;
	int error;
	int start_exon_error;
	char gene[100];
	genome_chunk translation;
	genome_chunk *utr5;
	int utr5_cnt;
	genome_chunk *utr3;
	int utr3_cnt;
	genome_chunk *exons;
	int exon_cnt;
	Strand strand;
	char biotype[512];
};

struct curl_data {
	char *memory;
 	size_t size;
};
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
struct vcf_mut_info {
        char mutation[512];
        int pos;
	int alt_pos;
        int offset;
        int flag;
};
struct contig_info {
        char id[50];
        int len;
        int mut_cnt;
        int *deleted_pos;
        int del_cnt;
        vcf_mut_info *mut;
};
struct transcriptome_mutations {
	contig_info *contigs;
	int contig_cnt;
};
struct exon {
        int g_st;
        int g_en;
        int len;
        int transcript_st;
        int transcript_en;
};
struct mut_info {
        int coordinate;
        int present;
	int indel;
        char from[100];
        char to[100];
};
struct genome_info {
	char transcript_id[100];
	char ref_id[100];
        char chr[50];
	char class_code;
        char gene[512];
	char biotype[512];
	int start;
	int end;
	char *embl_coords;
	char *alt_coords;
        Strand strand;
	char metadata[512];
	char category[512];
	char **muts;
	int mut_cnt;
};
struct transcriptome_info {
	int transcript_cnt;
	genome_info *transcripts;
};
struct seq_info {
	char seq[100];
	char no_ptm[100];
	transcriptome_info *transcriptomes;
	int uniprot_cnt;
        char **uniprot_accessions;	
};

int count_seq(FILE *);
seq_info *get_seqs(FILE *, int);
void remove_ptm(char *, char *);
void get_uniprot(FILE *, seq_info **, int);
int get_transcripts(FILE *, seq_info **, int, int, transcriptome_mutations *);
int sort_alpha(const void *, const void *);
void make_translation(char ***, char *);
void parse_title(genome_info *, char *, int, int, transcriptome_mutations *);
char *move_ptr(char *,char);
void check_mut(exon *, int, mut_info *, int);
int parse_exons(exon **, char *, char *, int);
int transcript_to_genome_loc(exon *, int, int, int);
void get_metadata(genome_info *, embl_info);
int count_alt_vcf(FILE *, int *);
vcf_info *read_alt_vcf(FILE *, contig_info **, int, int);
void create_mut_table(vcf_info *, int, contig_info *, int);
void add_deleted_pos_memory(int **, int, int);
static size_t receive(void *, size_t, size_t, void *);
int liftover(int, contig_info);
void get_embl(seq_info **, int, transcriptome_mutations *, int); 
void get_canonical_frame_info(embl_info *);
char *calculate_info(genome_info *, embl_info);
Region find_region(int, embl_info);
int find_frame(int, embl_info);
void get_coords_string(char **, embl_info, Strand);
tracking_info *get_tracking(FILE *, int);
void add_tracking_info(seq_info **,int,int,tracking_info *,int);
void add_mutation_info(genome_info *,vcf_mut_info **,int);
void make_mutation_strings(char ***, int, vcf_mut_info *, int);
void reset_flags(contig_info *);
vcf_mut_info *filter_mut_by_exons(exon *, int, contig_info *, int *);
void print_help();
/*
origins.c: try to work out where transcriptome-based peptides came from - totally new transcript? intron? 'UTR'?
Searches sequences in uniprot - if there, classed 'conventional' and added to _origins_prot.csv output
Search in translated transcriptome; get metadata from class_code and ensembl
Also calculates liftover for alternate transcriptomes and whether sequences include variant sites 
*/
/**********************************************************************************************************/
void main(int argc,char **argv)
{
	FILE *f, *g, *h, *e, *file_out_p, *file_out_t, *tracking_file, *discard_file, *unco_file;
	char output1[532], output2[532], stem[512], **vcf_files = NULL, all[500], **alt_tr_files = NULL, pep_file[512], prot_file[512], norm_tr_file[512], c = '\0', *pnt = NULL, *ptr = NULL, *loc = NULL, tracking_fn[512], discard_fn[532], unco_fn[532];
	int i, j, k, m, cnt = 0, num = 0, *transcript_cnt = NULL, vcf_lines_cnt = 0, contig_cnt = 0, vcf_cnt = 0, alt_cnt = 0, max, tracking_cnt = 0;
	seq_info *seqs = NULL;
	transcriptome_mutations *all_mutations = NULL;  
        vcf_info *vcf_lines = NULL;
	tracking_info *tracking = NULL;
	Mode_type mode = EMBL_ON;

	pep_file[0] = '\0';
	prot_file[0] = '\0';
	norm_tr_file[0] = '\0';
	discard_fn[0] = '\0';
	
	printf("Running origins with following options:\n");
	while ((c = getopt (argc, argv, "v:n:a:p:d:t:r:sh")) != -1) {
                switch (c) {
			case 'r':
				tracking_fn[0] = '\0';
				strcpy(tracking_fn, optarg);
				printf("-r\ttracking file %s\n", tracking_fn);
				break;
                        case 'v':
				all[0] = '\0';
				strcpy(all, optarg);
				ptr = all;
				while ((loc = strstr(ptr, ",")) != NULL) {
					vcf_cnt++;
					ptr = loc + 1;
				}
				vcf_cnt++;
                		if ((vcf_files=calloc(vcf_cnt, sizeof(char *)))==NULL) {
                        		printf("memory allocation error in main\n");
                        		exit(0);
                		}
				ptr = all;
								
				for (i = 0; (pnt = move_ptr(ptr,',')) != NULL; i++) {
                			if ((vcf_files[i]=calloc(strlen(ptr) + 1, sizeof(char)))==NULL) {
                        			printf("memory allocation error in main\n");
                        			exit(0);
                			}
					strcpy(vcf_files[i],ptr);
					printf("-v\tvcf file %s\n", ptr);
					*pnt = ',';
					ptr = pnt+1;
				}
                		if ((vcf_files[i]=calloc(strlen(ptr) + 1, sizeof(char)))==NULL) {
                        		printf("memory allocation error in main\n");
                        		exit(0);
                		}
				strcpy(vcf_files[i],ptr);
				printf("-v\tvcf file %s\n", ptr);
				break;
			case 'a': 
				all[0] = '\0';
				strcpy(all, optarg);
				ptr = all;
				while ((loc = strstr(ptr, ",")) != NULL) {
					alt_cnt++;
					ptr = loc + 1;
				}
				alt_cnt++;
                		if ((alt_tr_files=calloc(alt_cnt, sizeof(char *)))==NULL) {
                        		printf("memory allocation error in main\n");
                        		exit(0);
                		}
				ptr = all;
								
				for (i = 0; (pnt = move_ptr(ptr,',')) != NULL; i++) {
                			if ((alt_tr_files[i]=calloc(strlen(ptr) + 1, sizeof(char)))==NULL) {
                        			printf("memory allocation error in main\n");
                        			exit(0);
                			}
					strcpy(alt_tr_files[i],ptr);
					printf("-a\talt transcriptome file %s\n", ptr);
					*pnt = ',';
					ptr = pnt+1;
				}
                		if ((alt_tr_files[i]=calloc(strlen(ptr) + 1, sizeof(char)))==NULL) {
                        		printf("memory allocation error in main\n");
                        		exit(0);
                		}
				strcpy(alt_tr_files[i],ptr);
				printf("-a\talt transcriptome file %s\n", ptr);
				break;
			case 'p':
				strcpy(pep_file, optarg);
				printf("-p\tpeptide file %s\n", pep_file);
				break;
			case 'd':
				strcpy(prot_file, optarg);
				printf("-d\tprotein database %s\n", prot_file);
				break;
			case 'n':
				strcpy(norm_tr_file, optarg);
				printf("-n\t'normal' transcriptome file (no mutations)%s\n", norm_tr_file);
				break;
			case 's':
                               mode = SIMPLE;
                               printf("-s\tsimple mode chosen\n\t\t(Ensembl will not be consulted for detailed origins of sequences.)\n");
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
	fflush(stdout);
	if (vcf_cnt != alt_cnt) {
		printf("Provide a vcf file for each alt transcriptome file, in comma-separated lists.\n\te.g. -v v1.vcf,v2.vcf -a alt1.fa,alt2.fa\n\n");
		print_help();
		exit(0);
	}
	if (pep_file[0] == '\0' || prot_file[0] == '\0' || norm_tr_file[0] == '\0' || tracking_fn[0] == '\0') {
		printf("Specify .txt file of peptides with -p, .fa transcriptome (without mutations) file with -n, .fa uniprot db with -d, and .tracking file (from cuffcompare) with -r. If there are alt transcriptomes including mutations, specify vcf files in comma-separated lists with -v, and corresponding transcriptomes in comma-separated lists with -a\n\n");
		print_help();
		exit(0);
	}
        if ((transcript_cnt=calloc(vcf_cnt + 1, sizeof(int)))==NULL) {
               	printf("memory allocation error in main\n");
                exit(0);
        }
	if (vcf_cnt > 0) {
                if ((all_mutations=calloc(vcf_cnt, sizeof(transcriptome_mutations)))==NULL) {
                     	printf("memory allocation error in main\n");
                	exit(0);
                }
//create mutation tables for liftover (i.e. contigs array)
		for (i = 0; i < vcf_cnt; ++i) {
			if ((f = fopen(vcf_files[i],"r")) == NULL) {
                		printf("Can't open file %s\n",vcf_files[i]);
                		exit(0);
        		}
			all_mutations[i].contig_cnt = 0;
			vcf_lines_cnt = count_alt_vcf(f, &all_mutations[i].contig_cnt);
        		printf("Found %d contig lines and %d lines of mutation information in %s\n", all_mutations[i].contig_cnt, vcf_lines_cnt, vcf_files[i]);
        		vcf_lines = read_alt_vcf(f, &all_mutations[i].contigs, vcf_lines_cnt, all_mutations[i].contig_cnt);
			for (j = 0; j < all_mutations[i].contig_cnt; ++j) {
				if (all_mutations[i].contigs[j].mut_cnt > 0) {
                			if ((all_mutations[i].contigs[j].mut=calloc(all_mutations[i].contigs[j].mut_cnt, sizeof(vcf_mut_info)))==NULL) {
                	        		printf("memory allocation error in main\n");
                	        		exit(0);
                			}
				}
			}
        		create_mut_table(vcf_lines, vcf_lines_cnt, all_mutations[i].contigs, all_mutations[i].contig_cnt);
			if (vcf_lines != NULL)
                		free(vcf_lines);
			vcf_lines_cnt = 0;
			fclose(f);
		}
/*		j=1;
		for (i=0; i<all_mutations[j].contig_cnt; ++i) {
       			for (k=0; k<all_mutations[j].contigs[i].mut_cnt; k++)
        	        	printf("%d.%d) %s, %d, %d\n", i, k, all_mutations[j].contigs[i].mut[k].mutation, all_mutations[j].contigs[i].mut[k].pos, all_mutations[j].contigs[i].mut[k].offset);
        	}*/
	}
	if ((f = fopen(pep_file,"r")) == NULL) {
                printf("Can't open file %s\n",pep_file);
                exit(0);
        }
	if ((h = fopen(prot_file,"r")) == NULL) {
                printf("Can't open file %s\n",prot_file);
                exit(0);
        }
	strcpy(stem, pep_file);
	if ((pnt = strstr(stem, ".txt")) == NULL) {
		printf("Specify .txt file of peptides, .fa transcriptome file and .fa uniprot db\n");
		exit(0);
	}
	*pnt = '\0';
	//get output files ready
	sprintf(output1, "%s_origins_prot.csv", stem);
	sprintf(output2, "%s_origins_rna.csv", stem);
	sprintf(discard_fn, "%s_origins_discard.txt", stem);
	if ((file_out_p = fopen(output1,"w")) == NULL) {
                printf("Can't open file %s\n",output1);
                exit(0);
        }
	if ((file_out_t = fopen(output2,"w")) == NULL) {
                printf("Can't open file %s\n",output2);
                exit(0);
        }
	if ((discard_file = fopen(discard_fn,"w")) == NULL) {
                printf("Can't open file %s\n",discard_fn);
                exit(0);
        }
	
	cnt = count_seq(f);
	printf("Seq: %d\n", cnt);	
	fflush(stdout);
	seqs = get_seqs(f, cnt);
	fclose(f);

	qsort(seqs, cnt, sizeof(seq_info), sort_alpha);
	
	for (i = 0; i < cnt; ++i) {
		remove_ptm(seqs[i].no_ptm, seqs[i].seq);
	}
//get nearest_ref and class_code for transcripts
	if ((tracking_file = fopen(tracking_fn,"r")) == NULL) {
                printf("Can't open file %s\n",tracking_fn);
                exit(0);
        }
	tracking_cnt = count_seq(tracking_file);
	tracking = get_tracking(tracking_file, tracking_cnt);

	fclose(tracking_file);
//get uniprot matches	
	get_uniprot(h, &seqs, cnt);
	fclose(h);
//get transcript matches
	for (i = 0; i < cnt; ++i) {
		if ((seqs[i].transcriptomes=calloc(vcf_cnt + 1, sizeof(transcriptome_info)))==NULL) {
             		printf("memory allocation error in main\n");
               		exit(0);
        	}
	}
	for (i = 0; i < vcf_cnt + 1; i++) {
		if (i == 0) {
			if ((e = fopen(norm_tr_file,"r")) == NULL) {
                		printf("Can't open file %s\n",norm_tr_file);
                		exit(0);
        		}
			printf("searching transcripts for matches with peptide sequences in %s...\n", norm_tr_file);
			fflush(stdout);
			transcript_cnt[i] = get_transcripts(e, &seqs, cnt, i, NULL);
			printf("done - %d transcripts searched\n", transcript_cnt[i]);
			fflush(stdout);
		} else {
			if ((e = fopen(alt_tr_files[i-1],"r")) == NULL) {
                		printf("Can't open file %s\n",alt_tr_files[i-1]);
                		exit(0);
        		}
			printf("searching transcripts for matches with peptide sequences in %s...\n", alt_tr_files[i-1]);
			fflush(stdout);
			transcript_cnt[i] = get_transcripts(e, &seqs, cnt, i, &(all_mutations[i-1]));
			printf("done - %d transcripts searched\n", transcript_cnt[i]);
			fflush(stdout);
		}
		fclose(e);
	}
	printf("Transcriptome sequences searched:\n");
	for (i = 0; i < vcf_cnt + 1; i++)
		if (i == 0)
			printf("%d) %d from %s\n", i, transcript_cnt[i], norm_tr_file);
		else
			printf("%d) %d from %s\n", i, transcript_cnt[i], alt_tr_files[i-1]);
	fflush(stdout);
	for (i = 0; i < vcf_cnt; ++i) {
		if (vcf_files[i] != NULL)
			free(vcf_files[i]);
	}
	if (vcf_files != NULL)
		free(vcf_files);
	if (transcript_cnt != NULL)
		free(transcript_cnt);

//Add tracking info
	add_tracking_info(&seqs,cnt,vcf_cnt +1,tracking,tracking_cnt);	
	if (tracking != NULL)
		free(tracking);
//Add mutation info
	if (mode == EMBL_ON) {
               get_embl(&seqs, cnt, all_mutations, vcf_cnt+1);
       	}

	for (i = 0; i < vcf_cnt; ++i) {
		for (j = 0; j < all_mutations[i].contig_cnt; ++j) {
			if (all_mutations[i].contigs[j].mut != NULL) 
				free(all_mutations[i].contigs[j].mut);
			if (all_mutations[i].contigs[j].deleted_pos != NULL)
				free(all_mutations[i].contigs[j].deleted_pos);
		}
		if (all_mutations[i].contigs != NULL)
			free(all_mutations[i].contigs);
	}
	if (all_mutations != NULL)
		free(all_mutations);
	
	fprintf(file_out_p, "Num,Sequence,\"Protein Accession\"\n");
	for (i = 0; i < cnt; ++i) {
		for (j = 0; j < seqs[i].uniprot_cnt; ++j) {
			fprintf(file_out_p, "%d,%s,\"%s\"\n", i+1,seqs[i].seq,seqs[i].uniprot_accessions[j]);
		}
	}
	fclose(file_out_p); //output prot

	//in SIMPLE mode, print out unconventional peptide list
	if (mode == SIMPLE) {
		sprintf(unco_fn, "%s_origins_unconventional.txt", stem);
		if ((unco_file = fopen(unco_fn,"w")) == NULL) {
                	printf("Can't open file %s\n",unco_fn);
                	exit(0);
        	}
		for (i = 0; i < cnt; ++i) {
			max = 0;
                	for (j = 0; j < vcf_cnt + 1; ++j) {
                        	if (seqs[i].transcriptomes[j].transcript_cnt > max)
                                max = seqs[i].transcriptomes[j].transcript_cnt;
                	}
			if (!seqs[i].uniprot_cnt && max)
				fprintf(unco_file, "%s\n",seqs[i].seq);
		}
		fclose(unco_file);
	}

	fprintf(file_out_t, ",,,");
	for (i = 0; i < vcf_cnt+1; i++)
		fprintf(file_out_t, "%s,,,,,,,,,,,,", (i==0) ? norm_tr_file : alt_tr_files[i-1]);
	fprintf(file_out_t, "\n");
	fprintf(file_out_t, "Num,Sequence,Conventional explanation?,");
	for (i = 0; i < vcf_cnt+1; i++)	{
		fprintf(file_out_t, "Transcript,Mutations,Metadata,Category,Chromosome,\"Peptide coordinates (from transcript)\",Strand,\"Reference id\",\"Alt genome coordinates\",\"Ensembl coordinates\",Biotype,Gene");
		if (i < vcf_cnt)
			fprintf(file_out_t,",");
	}
	fprintf(file_out_t, "\n");
	
	for (i = 0; i < cnt; ++i) {
		max = 0;
		for (j = 0; j < vcf_cnt + 1; ++j) {
			if (seqs[i].transcriptomes[j].transcript_cnt > max)
				max = seqs[i].transcriptomes[j].transcript_cnt;
		}
		if (!max && !seqs[i].uniprot_cnt) { 
			fprintf(discard_file, "%s\n",seqs[i].seq);
		}
		for (k = 0; k < max; ++k) {
			fprintf(file_out_t, "%d,%s,%s,", i+1,seqs[i].seq,(seqs[i].uniprot_cnt) ? "Y" : "N");
			for (j = 0; j < vcf_cnt + 1; ++j) {
				if (k < seqs[i].transcriptomes[j].transcript_cnt) {
					fprintf(file_out_t, "%s,", seqs[i].transcriptomes[j].transcripts[k].transcript_id);
					if (seqs[i].transcriptomes[j].transcripts[k].mut_cnt) {
						for (m = 0; m < seqs[i].transcriptomes[j].transcripts[k].mut_cnt; ++m) {
							fprintf(file_out_t, "%s",seqs[i].transcriptomes[j].transcripts[k].muts[m]);
							if (m < seqs[i].transcriptomes[j].transcripts[k].mut_cnt - 1)
								fprintf(file_out_t, ";");
						}
					}
					fprintf(file_out_t, ",\"%s\",%s,%s,%d-%d,", seqs[i].transcriptomes[j].transcripts[k].metadata,seqs[i].transcriptomes[j].transcripts[k].category,seqs[i].transcriptomes[j].transcripts[k].chr, seqs[i].transcriptomes[j].transcripts[k].start, seqs[i].transcriptomes[j].transcripts[k].end);
					if (seqs[i].transcriptomes[j].transcripts[k].strand == PLUS) {
						fprintf(file_out_t, "PLUS,");
					} else if (seqs[i].transcriptomes[j].transcripts[k].strand == MINUS) {
						fprintf(file_out_t, "MINUS,");
					} else {
						fprintf(file_out_t, "UNKNOWN,");
					}
					if (seqs[i].transcriptomes[j].transcripts[k].class_code != 'u') {
                               			fprintf(file_out_t, "%s,", seqs[i].transcriptomes[j].transcripts[k].ref_id);
						if (seqs[i].transcriptomes[j].transcripts[k].alt_coords != NULL)
                               				fprintf(file_out_t, "%s", seqs[i].transcriptomes[j].transcripts[k].alt_coords);
						fprintf(file_out_t, ",");
						if (seqs[i].transcriptomes[j].transcripts[k].embl_coords != NULL)
                               				fprintf(file_out_t, "%s", seqs[i].transcriptomes[j].transcripts[k].embl_coords);
                               			fprintf(file_out_t, ",\"%s\",%s,", seqs[i].transcriptomes[j].transcripts[k].biotype, seqs[i].transcriptomes[j].transcripts[k].gene);
					} else {
						fprintf(file_out_t, ",,,,,");
					}
				} else {
					fprintf(file_out_t, ",,,,,,,,,,,,");
				}
			}
			fprintf(file_out_t,"\n");
		}
	}
//free and cleanup
	for (i = 0; i < alt_cnt; ++i) {
		if (alt_tr_files[i] != NULL)
			free(alt_tr_files[i]);
	}
	if (alt_tr_files != NULL)
		free(alt_tr_files);

	for (i = 0; i < cnt; ++i) {
		if (seqs[i].uniprot_cnt) {
			for (j = 0; j < seqs[i].uniprot_cnt; ++j) {
				if (seqs[i].uniprot_accessions[j] != NULL)
					free(seqs[i].uniprot_accessions[j]);
			}
			if (seqs[i].uniprot_accessions != NULL)
				free(seqs[i].uniprot_accessions);
		}
		for (j=0; j < vcf_cnt + 1; j++) {
			for (k = 0; k < seqs[i].transcriptomes[j].transcript_cnt; ++k) {
				if (seqs[i].transcriptomes[j].transcripts[k].embl_coords != NULL)
					free(seqs[i].transcriptomes[j].transcripts[k].embl_coords);
				if (seqs[i].transcriptomes[j].transcripts[k].alt_coords != NULL)
					free(seqs[i].transcriptomes[j].transcripts[k].alt_coords);
				for (m=0; m < seqs[i].transcriptomes[j].transcripts[k].mut_cnt; ++m) {
					if (seqs[i].transcriptomes[j].transcripts[k].muts[m] != NULL)
						free(seqs[i].transcriptomes[j].transcripts[k].muts[m]);
				}
				if (seqs[i].transcriptomes[j].transcripts[k].muts != NULL)
					free(seqs[i].transcriptomes[j].transcripts[k].muts);
			}
			if (seqs[i].transcriptomes[j].transcripts != NULL)
				free(seqs[i].transcriptomes[j].transcripts);
		}
		if (seqs[i].transcriptomes != NULL)
			free(seqs[i].transcriptomes);
	}
	if (seqs != NULL)
		free(seqs);	

	fclose(file_out_t); //output transcriptome
	fclose(discard_file); //file of peptides to discard because no source in transcriptomes or protein db (presumably artificial junction peptides)
	printf("Finished\n");	
	exit(0);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void print_help() {
	printf("Required Input:\n\t-p\ttxt file listing peptide sequences to check\n");
	printf("\t-r\tCuffcompare output tracking file\n"); 
	printf("\t-d\tstandard protein database (e.g. Uniprot)\n");
	printf("\t-n\t‘normal’ transcriptome i.e. without variants added\n\t\t(nucleotide sequences in fasta format)\n");
	printf("Optional*:\n");
	printf("\t-a\tcomma-separated list of ‘alternate’ transcriptomes\n\t\t(i.e. without variants added (nucleotide sequences in fasta format)\n");
	printf("\t-v\tcomma-separated list of the vcf files used to produce the corresponding transcriptomes\n");
	printf("\t-s\tchoose simple mode - no detailed info directly from Ensembl\n\t\tThis may be useful for a quick look at results or if the Ensembl REST API is malfunctioning)\n");
	printf("\t-h\tprint help\n* -a and –v are optional, but if used the files should correspond;\n\t(e.g. –a transcriptome1.fa,transcriptome2.fa –v v1.vcf,v2.vcf)\n");
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int count_seq(FILE *f)
{
	int cnt = 0;
	char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (isalpha(line[0]))    
                    cnt++;
        }
        rewind(f);

        return(cnt);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
seq_info *get_seqs(FILE *f, int cnt)
{
	seq_info *tmp = NULL;
	char line[MAXSTR];
	int i = 0;

	if ((tmp = calloc(cnt, sizeof(seq_info))) == NULL) {
		printf("Memory allocation error in get_seqs()\n");
                exit(0);
	}

        while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (isalpha(line[0])) {  
			while (line[strlen(line) -1] == '\n' || line[strlen(line) -1] == '\r' || line[strlen(line) -1] == ' ')
				line[strlen(line) -1] = '\0';
			strcpy(tmp[i].seq, line);
			i++; 
		}
        }
	return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
tracking_info *get_tracking(FILE *f, int cnt)
{
	int i=0;
        char line[MAXSTR], ref[512], *ptr = NULL, *pnt = NULL;
        tracking_info *tmp;

        if ((tmp = calloc(cnt,sizeof(tracking_info))) == NULL) {
                printf("Memory allocation error in get_tracking()\n");
                exit(0);
        }

        while(fgets(line,MAXSTR - 1,f) != NULL) {
           	ref[0] = '\0';
                sscanf(line, "%s %*[^\t] %s %c", tmp[i].cuff_transcript, ref, &tmp[i].class_code);
                if (ref[0] != '-') {
			ptr = ref;
			pnt = move_ptr(ref,'|');
			if (pnt == NULL) {
				printf("Unexpected formatting of nearest_ref in tracking file: %s\n", ref);
				exit(0);
			}	
			ptr = pnt + 1;
			strcpy(tmp[i].nearest_ref, ptr);
                } else {
			strcpy(tmp[i].nearest_ref, ref);
		}
                i++;
        }

        if (i != cnt) {
		printf("Numbers don't match in get_tracking - %d vs %d", i, cnt);
		exit(0);
        }
        return tmp;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void add_tracking_info(seq_info **seqs,int cnt,int transcriptome_cnt,tracking_info *tracking,int tracking_cnt)
{
	int i,j,k,m;

	for (i = 0; i < cnt; ++i) {
		for (j=0; j < transcriptome_cnt; ++j) {
			for (k=0; k < (*seqs)[i].transcriptomes[j].transcript_cnt; ++k) {
				for (m = 0; m < tracking_cnt; ++m) {
                                	if (strcmp(tracking[m].cuff_transcript, (*seqs)[i].transcriptomes[j].transcripts[k].transcript_id) == 0) {
						(*seqs)[i].transcriptomes[j].transcripts[k].class_code = tracking[m].class_code;
						strcpy((*seqs)[i].transcriptomes[j].transcripts[k].ref_id, tracking[m].nearest_ref);
						break;
					}

				}
			}
		}
	}
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void add_mutation_info(genome_info *transcript,vcf_mut_info **filtered_mut,int filtered_cnt)
{
	int i, tmp_mut_cnt = 0;

	//get mutations for peptide
	for (i = 0; i < filtered_cnt; ++i) {
		if (transcript->strand == PLUS || transcript->strand == UNKNOWN) {
			if ((transcript->start -1 < (*filtered_mut)[i].alt_pos) && (transcript->end +1 > (*filtered_mut)[i].alt_pos)) {
				tmp_mut_cnt++;
				(*filtered_mut)[i].flag = 1;
			}
		} else {
			if ((transcript->end -1 < (*filtered_mut)[i].alt_pos) && (transcript->start +1 > (*filtered_mut)[i].alt_pos)) {
				tmp_mut_cnt++;
				(*filtered_mut)[i].flag = 1;
			}	
		}
	}
	//into the transcript, save the mutation info we just found
	if (tmp_mut_cnt) {
		transcript->mut_cnt = tmp_mut_cnt;
		make_mutation_strings(&(transcript->muts), tmp_mut_cnt, *filtered_mut, filtered_cnt);
	}
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void make_mutation_strings(char ***muts, int cnt, vcf_mut_info *filtered_mut, int filtered_cnt)
{
	char **tmp = NULL;
	int i,j = 0, sanity_check = 0;

	if ((tmp = calloc(cnt, sizeof(char *))) == NULL) {
         	printf("memory allocation error in make_mutation_strings()\n");
        	exit(0);
	}
	for (i = 0; i < filtered_cnt; ++i) {
		if (filtered_mut[i].flag) {
			sanity_check++;
			if (sanity_check > cnt) {
				printf("Counting problem in make_mutation_strings()\n");
				exit(0);
			}
			if ((tmp[j] = calloc(strlen(filtered_mut[i].mutation) + 1, sizeof(char))) == NULL) {
         			printf("memory allocation error in make_mutation_strings()\n");
        			exit(0);
			}
			strcpy(tmp[j],filtered_mut[i].mutation);
			j++;
		}	
	}
	*muts = tmp;
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void reset_flags(contig_info *c)
{
        int i;
        for (i=0; i < c->mut_cnt; ++i)
                c->mut[i].flag = 0;
        return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void remove_ptm(char *no_ptm, char *pep)
{
        char *pnt = NULL, *ptr = NULL;
	
	while (*pep != '\0') {
		if (*pep == '(') {
			while (*pep != ')')
				pep++;
			pep++;
		}
		else
			*no_ptm++ = *pep++;
	}
	*no_ptm = '\0'; 
        return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void get_canonical_frame_info(embl_info *embl)
{
	int i, j, start_exon;

//for each exon, check if translation start fits in it. index = start_exon.
//calculate bases_from_start for each exon as appropriate
	if (embl->translation.st) {
		start_exon = -1;
		if (embl->strand == PLUS) {
			for (j = 0; j < embl->exon_cnt; j++) {
				if ((embl->translation.st < embl->exons[j].en + 1) && (embl->translation.st > embl->exons[j].st - 1))
					start_exon = j;
			}
			if (start_exon < 0) {
				printf("start_exon not found - PLUS strand\n");
				embl->start_exon_error = 1;
			} else {
				if (start_exon < (embl->exon_cnt - 1))
					embl->exons[start_exon + 1].bases_from_start = embl->exons[start_exon].en + 1 - embl->translation.st;
				for (j = start_exon + 2; j < embl->exon_cnt; j++) 
					embl->exons[j].bases_from_start = embl->exons[j - 1].bases_from_start + embl->exons[j - 1].en + 1 - embl->exons[j - 1].st;
				embl->exons[start_exon].bases_from_start = embl->exons[start_exon].st - embl->translation.st;
				for (j = start_exon - 1; j > -1; j--) 
					embl->exons[j].bases_from_start = embl->exons[j + 1].bases_from_start - (embl->exons[j].en + 1 - embl->exons[j].st);
			}
		} else if (embl->strand == MINUS) {
			for (j = 0; j < embl->exon_cnt; j++) {
				if ((embl->translation.en < embl->exons[j].en + 1) && embl->translation.en > embl->exons[j].st - 1)
					start_exon = j;
			}
			if (start_exon < 0) {
				printf("start_exon not found - MINUS strand\n");
				embl->start_exon_error = 1;
			} else {
				if (start_exon < (embl->exon_cnt - 1))
					embl->exons[start_exon + 1].bases_from_start = embl->translation.en - embl->exons[start_exon].st + 1;
				for (j = start_exon + 2; j < embl->exon_cnt; j++) 
					embl->exons[j].bases_from_start = embl->exons[j-1].bases_from_start + (embl->exons[j-1].en - embl->exons[j-1].st + 1);	
				embl->exons[start_exon].bases_from_start = embl->translation.en - embl->exons[start_exon].en;	
				for (j = start_exon - 1; j > -1; j--) 
					embl->exons[j].bases_from_start = embl->exons[j].st - embl->exons[j].en - 1 + embl->exons[j+1].bases_from_start;	
			}
		} else {
			printf("Found 'translated' unknown strand in embl. Freak out!\n");
			exit(0);
		}
	}
	
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void get_embl(seq_info **seqs, int cnt, transcriptome_mutations *all_mutations, int transcriptome_cnt) 
{
	embl_info *embls = NULL, *realloc_tmp = NULL;
	CURL *handle;
        CURLcode result;
        curl_data data;
        cJSON *json = NULL, *object_type = NULL, *translation = NULL, *parent = NULL, *utr = NULL, *utrs = NULL, *exons = NULL, *exon = NULL, *start = NULL, *end = NULL, *strand = NULL, *biotype = NULL, *error = NULL, *version = NULL;
	int i, j, k, m, index = -1, match = -1, embl_cnt = 0, error_cnt = 0, good = 0;
	char url[512], *string, trimmed[50], *pnt = NULL;
	genome_chunk *realloc_chunk = NULL;

	curl_global_init(CURL_GLOBAL_ALL);
	handle = curl_easy_init();
	curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, receive);
        curl_easy_setopt(handle, CURLOPT_WRITEDATA, (void *)&data);

	printf("Getting ensembl info for reference transcripts...\n");
	fflush(stdout);

//set up first embl_info
	if ((embls = calloc(1, sizeof(embl_info))) == NULL) {
         	printf("memory allocation error in get_embl()\n");
        	exit(0);
	}
						
	for (i = 0; i < cnt; ++i) {
		for (j = 0; j < transcriptome_cnt; ++j) {
			for (k = 0; k < (*seqs)[i].transcriptomes[j].transcript_cnt; ++k) {
				index = -1; 
				if ((*seqs)[i].transcriptomes[j].transcripts[k].class_code != 'u') {
					trimmed[0] = '\0';
					strcpy(trimmed, (*seqs)[i].transcriptomes[j].transcripts[k].ref_id);
					if ((pnt = strchr(trimmed, '.')) != NULL)
						*pnt = '\0';
					for (m = 0; m < embl_cnt; m++) {
						if ((strcmp(trimmed, embls[m].id)) == 0) 
							index = m;
					}
					if (index < 0) { //get new info from ensembl about this transcript
						index = embl_cnt;
	//set up embl_info space
						if (embl_cnt) {
                                        		if ((realloc_tmp = realloc(embls, ((embl_cnt + 1) * sizeof(embl_info)))) == NULL) {
                                        	               	printf("memory allocation error in get_embl()\n");
                                        	               	exit(0);
                                        	        }
							embls = realloc_tmp;
							memset(embls + embl_cnt, 0, sizeof(embl_info));	
						}
	//add info
	//id		
						strcpy(embls[embl_cnt].id, trimmed);
						embls[embl_cnt].id_num = (int)atoi(embls[embl_cnt].id +4);
						data.memory = malloc(1);
        					data.size = 0;
						sprintf(url, "http://rest.ensembl.org/lookup/id/%s?expand=1;utr=1;content-type=application/json", embls[embl_cnt].id);
						curl_easy_setopt(handle, CURLOPT_URL, url);
						if ((result = curl_easy_perform(handle)) != CURLE_OK) {
                					printf("curl_easy_perform() failed: %s\n", curl_easy_strerror(result));
                					exit(0);
        					}
        					json = cJSON_Parse(data.memory);
						error = cJSON_GetObjectItemCaseSensitive(json, "error");
						if (cJSON_IsString(error) && (error->valuestring != NULL)) {
							embls[embl_cnt].error = 1;
							error_cnt++;
							printf("Failed to find info on %s\n", embls[embl_cnt].id);
						} else {
							good++;
	//gene					
							parent = cJSON_GetObjectItemCaseSensitive(json, "Parent");
    							if (cJSON_IsString(parent) && (parent->valuestring != NULL)) {
								strcpy(embls[embl_cnt].gene, parent->valuestring);
							} else {
                						printf("getting gene name failed for transcript %s\n", embls[embl_cnt].id);
                						exit(0);
							}
	//version
							version = cJSON_GetObjectItemCaseSensitive(json, "version");
    							if (cJSON_IsNumber(version) && (version->valueint > -1)) {
								embls[embl_cnt].version = version->valueint;
							} else {
                						printf("no version number for transcript %s\n", embls[embl_cnt].id);
                						exit(0);
							}
	//translation (if present)
							translation = cJSON_GetObjectItemCaseSensitive(json, "Translation");
        						if (cJSON_IsObject(translation)) {
                						start = cJSON_GetObjectItemCaseSensitive(translation, "start");
                						end = cJSON_GetObjectItemCaseSensitive(translation, "end");
                						if (cJSON_IsNumber(start) && cJSON_IsNumber(end)) {
									embls[embl_cnt].translation.embl_st = start->valueint;
									embls[embl_cnt].translation.embl_en = end->valueint;
								}
                					}
	//utrs (if present)	
							utrs = cJSON_GetObjectItemCaseSensitive(json, "UTR");
							cJSON_ArrayForEach(utr, utrs) {
                						start = cJSON_GetObjectItemCaseSensitive(utr, "start");
                						end = cJSON_GetObjectItemCaseSensitive(utr, "end");
								object_type = cJSON_GetObjectItemCaseSensitive(utr, "object_type");
								if ((strcmp(object_type->valuestring, "five_prime_UTR")) == 0) {
									if (!embls[embl_cnt].utr5_cnt) {
										if ((embls[embl_cnt].utr5 = calloc(1, sizeof(genome_chunk))) == NULL) {
              										printf("memory allocation error in get_embl()\n");
        										exit(0);
										}
									} else {
                                        					if ((realloc_chunk = realloc(embls[embl_cnt].utr5, (embls[embl_cnt].utr5_cnt + 1) * sizeof(genome_chunk))) == NULL) {
                                        	        	       			printf("memory allocation error in get_embl()\n");
                                        	        	       			exit(0);
                                        	        			}
										embls[embl_cnt].utr5 = realloc_chunk;
										memset(embls[embl_cnt].utr5 + embls[embl_cnt].utr5_cnt, 0, sizeof(genome_chunk));	
									}
                							if (cJSON_IsNumber(start) && cJSON_IsNumber(end)) {
										embls[embl_cnt].utr5[embls[embl_cnt].utr5_cnt].embl_st = start->valueint;
										embls[embl_cnt].utr5[embls[embl_cnt].utr5_cnt].embl_en = end->valueint;
									}
									embls[embl_cnt].utr5_cnt++;
								}
								else if ((strcmp(object_type->valuestring, "three_prime_UTR")) == 0) {
                                        	        	        if (!embls[embl_cnt].utr3_cnt) {
                                        	        	                if ((embls[embl_cnt].utr3 = calloc(1, sizeof(genome_chunk))) == NULL) {
                                        	        	                        printf("memory allocation error in get_embl()\n");
                                        	        	                        exit(0);
                                        	        	                }
                                        	        	        } else {
                                        	                        	if ((realloc_chunk = realloc(embls[embl_cnt].utr3, (embls[embl_cnt].utr3_cnt + 1) * sizeof(genome_chunk))) == NULL) {
                                        	                        	        printf("memory allocation error in get_embl()\n");
                                        	                        	        exit(0);
                                        	                        	}
                                        	                        	embls[embl_cnt].utr3 = realloc_chunk;
                                        	                        	memset(embls[embl_cnt].utr3 + embls[embl_cnt].utr3_cnt, 0, sizeof(genome_chunk));
                                        	                	}
                                        	                	if (cJSON_IsNumber(start) && cJSON_IsNumber(end)) {
										embls[embl_cnt].utr3[embls[embl_cnt].utr3_cnt].embl_st = start->valueint;
										embls[embl_cnt].utr3[embls[embl_cnt].utr3_cnt].embl_en = end->valueint;
                                        	                	}
									embls[embl_cnt].utr3_cnt++;
								} else {
									printf("What type of utr is this?\n");
									exit(0);
								}
							}
	//strand					
							strand = cJSON_GetObjectItemCaseSensitive(json, "strand");
							if (cJSON_IsNumber(strand)) {
								if (strand->valueint == 1)
									embls[embl_cnt].strand = PLUS;
								else if (strand->valueint == -1)
									embls[embl_cnt].strand = MINUS;
								else {
									embls[embl_cnt].strand = UNKNOWN;
									printf("Strand unknown for %s", embls[embl_cnt].id);
								}
							}
	//biotype
							biotype = cJSON_GetObjectItemCaseSensitive(json, "biotype");
    							if (cJSON_IsString(biotype) && (biotype->valuestring != NULL)) {
								strcpy(embls[embl_cnt].biotype, biotype->valuestring);
							} else {
                						printf("getting biotype failed for transcript %s\n", embls[embl_cnt].id);
                						exit(0);
							}
	//exons	
							exons = cJSON_GetObjectItemCaseSensitive(json, "Exon");
							if (cJSON_IsArray(exons)) {		
								cJSON_ArrayForEach(exon, exons) {
                							start = cJSON_GetObjectItemCaseSensitive(exon, "start");
                							end = cJSON_GetObjectItemCaseSensitive(exon, "end");
									if (!embls[embl_cnt].exon_cnt) {
										if ((embls[embl_cnt].exons = calloc(1, sizeof(genome_chunk))) == NULL) {
              										printf("memory allocation error in get_embl()\n");
        										exit(0);
										}
									} else {
                	                        				if ((realloc_chunk = realloc(embls[embl_cnt].exons, (embls[embl_cnt].exon_cnt + 1) * sizeof(genome_chunk))) == NULL) {
                	                                       				printf("memory allocation error in get_embl()\n");
                	                                       				exit(0);
                	                                			}
										embls[embl_cnt].exons = realloc_chunk;
										memset(embls[embl_cnt].exons + embls[embl_cnt].exon_cnt, 0, sizeof(genome_chunk));	
									}
                							if (cJSON_IsNumber(start) && cJSON_IsNumber(end)) {
										embls[embl_cnt].exons[embls[embl_cnt].exon_cnt].embl_st = start->valueint;
										embls[embl_cnt].exons[embls[embl_cnt].exon_cnt].embl_en = end->valueint;
									}
									embls[embl_cnt].exon_cnt++;
								}
							} else {
								if ((embls[embl_cnt].exons = calloc(1, sizeof(genome_chunk))) == NULL) {
             								printf("memory allocation error in get_embl()\n");
        								exit(0);
								}
                						start = cJSON_GetObjectItemCaseSensitive(exons, "start");
                						end = cJSON_GetObjectItemCaseSensitive(exons, "end");
                						if (cJSON_IsNumber(start) && cJSON_IsNumber(end)) {
									embls[embl_cnt].exons[0].embl_st = start->valueint;
									embls[embl_cnt].exons[0].embl_en = end->valueint;
								}
								embls[embl_cnt].exon_cnt = 1;
							}
						}
						embl_cnt++;
						if (embl_cnt % 100 == 0) {
							printf("Info received for %d transcripts\n", embl_cnt); 
							fflush(stdout);
						}
        					cJSON_Delete(json);
						if (data.memory != NULL)
                					free(data.memory);
					}
//fill in info from embl in seqs transcript	
					if (!embls[index].error) {	
    						if (embls[index].gene[0] != '\0') {
							strcpy((*seqs)[i].transcriptomes[j].transcripts[k].gene,embls[index].gene);
						}
    						if (embls[index].biotype[0] != '\0') {
							strcpy((*seqs)[i].transcriptomes[j].transcripts[k].biotype,embls[index].biotype);
						}
	//to get 'normal' coordinates ensure en and st match embl_en and embl_st
						for (m = 0; m<embls[index].exon_cnt; ++m) {
							embls[index].exons[m].st = embls[index].exons[m].embl_st;
							embls[index].exons[m].en = embls[index].exons[m].embl_en;	
						}
						for (m = 0; m<embls[index].utr5_cnt; ++m) {
							 embls[index].utr5[m].st = embls[index].utr5[m].embl_st;
							 embls[index].utr5[m].en = embls[index].utr5[m].embl_en;
						}
						for (m = 0; m<embls[index].utr3_cnt; ++m) {
							embls[index].utr3[m].st = embls[index].utr3[m].embl_st;
							embls[index].utr3[m].en = embls[index].utr3[m].embl_en;
						}
						embls[index].translation.st = embls[index].translation.embl_st;
						embls[index].translation.en = embls[index].translation.embl_en;
	//get coords					
						get_coords_string(&((*seqs)[i].transcriptomes[j].transcripts[k].embl_coords), embls[index],(*seqs)[i].transcriptomes[j].transcripts[k].strand);
	//do liftover for proper en and st values if not 'normal' transcriptome
						if (j > 0) {
					//get right chromosome contig for right transcriptome:
							match = -1;
							for (m = 0; m < all_mutations[j-1].contig_cnt; ++m) {
                        					if (strcmp(all_mutations[j-1].contigs[m].id, (*seqs)[i].transcriptomes[j].transcripts[k].chr) == 0) {
                        						match = m;
								}
                					}
							if (match < 0) {
								printf("chromosome names don't match: alt transcriptome %d, chr %s\n", j-1, (*seqs)[i].transcriptomes[j].transcripts[k].chr);
								exit(0);
							}
					//do liftover
							for (m = 0; m < embls[index].exon_cnt; ++m) {
								embls[index].exons[m].st = liftover(embls[index].exons[m].embl_st, all_mutations[j-1].contigs[match]);
								embls[index].exons[m].en = liftover(embls[index].exons[m].embl_en, all_mutations[j-1].contigs[match]);
							}
							for (m = 0; m < embls[index].utr5_cnt; ++m) {
								embls[index].utr5[m].st = liftover(embls[index].utr5[m].embl_st, all_mutations[j-1].contigs[match]);
								embls[index].utr5[m].en = liftover(embls[index].utr5[m].embl_en, all_mutations[j-1].contigs[match]);
							}
							for (m = 0; m < embls[index].utr3_cnt; ++m) {
								embls[index].utr3[m].st = liftover(embls[index].utr3[m].embl_st, all_mutations[j-1].contigs[match]);
								embls[index].utr3[m].en = liftover(embls[index].utr3[m].embl_en, all_mutations[j-1].contigs[match]);
							}
							if (embls[index].translation.embl_st) {
								embls[index].translation.st = liftover(embls[index].translation.embl_st, all_mutations[j-1].contigs[match]);
								embls[index].translation.en = liftover(embls[index].translation.embl_en, all_mutations[j-1].contigs[match]);
							}
							get_coords_string(&((*seqs)[i].transcriptomes[j].transcripts[k].alt_coords), embls[index],(*seqs)[i].transcriptomes[j].transcripts[k].strand);
						}
	//fill in bases_from_start
						get_canonical_frame_info(&embls[index]);
						if (embls[index].start_exon_error) { 
							printf("\tseq transcript id %s, transcriptome %d\n",(*seqs)[i].transcriptomes[j].transcripts[k].transcript_id, j);
							printf("\tembls id %s, class_code %c\n", embls[index].id, (*seqs)[i].transcriptomes[j].transcripts[k].class_code);
							printf("\ttrans st %d\n", embls[index].translation.st);
							for (m = 0; m < embls[index].exon_cnt; ++m) {
								printf("\texon %d) st %d en %d\n", m+1, embls[index].exons[m].st,  embls[index].exons[m].en);
							}
						}
					}
				}
	//get metadata (for all seq, even if class_code == 'u')
				if ((*seqs)[i].transcriptomes[j].transcripts[k].class_code == 'u') {
					get_metadata(&((*seqs)[i].transcriptomes[j].transcripts[k]), embls[0]);
				} else {
					get_metadata(&((*seqs)[i].transcriptomes[j].transcripts[k]), embls[index]);
				}
			}
		}
	}
	
	printf("Downloaded info on %d transcripts from ensembl\nFailed to find info on %d %s\nProcessing...\n", good, error_cnt, (error_cnt == 1) ? "transcript" : "transcripts");
	fflush(stdout);

//free embls
	for (i = 0; i < embl_cnt; ++i) {
		if (!embls[i].error) {
			if (embls[i].exon_cnt && embls[i].exons != NULL)
				free(embls[i].exons);	
			if (embls[i].utr5_cnt && embls[i].utr5 != NULL)
				free(embls[i].utr5);	
			if (embls[i].utr3_cnt && embls[i].utr3 != NULL)
				free(embls[i].utr3);
		}
	}
	if (embls != NULL)
		free(embls);

	curl_easy_cleanup(handle);
	curl_global_cleanup();
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void get_coords_string(char **string, embl_info embl, Strand strand)
{
	char tmp[510], *ptr = NULL;
	
	tmp[0] = '\0';

	if (strand == PLUS || strand == UNKNOWN) {
		sprintf(tmp, "%d-%d", embl.exons[0].st, embl.exons[embl.exon_cnt-1].en);
	}	
	else { //== MINUS
		sprintf(tmp, "%d-%d", embl.exons[0].en, embl.exons[embl.exon_cnt-1].st);
	}
	if ((ptr = calloc(strlen(tmp) + 1, sizeof(char))) == NULL) {
        	printf("memory allocation error in get_coords_string()\n");
        	exit(0);
	}
	strcpy(ptr,tmp);
	*string = ptr;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void get_metadata(genome_info *transcript, embl_info embl)
{
	char *info, url[512], *string, trimmed[50], *pnt = NULL;

	switch(transcript->class_code) {
		case 'u':
			if ((info = calloc(22, sizeof(char))) == NULL) {
        			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			strcpy(info, "Utterly unknown/novel");
			strcpy(transcript->category, "novel");
			break;
		case 'x':
			if ((info = calloc(40 + strlen(transcript->ref_id), sizeof(char))) == NULL) {
              			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			sprintf(info, "exonic overlap with %s on opposite strand", transcript->ref_id);
			strcpy(transcript->category, "non-canonical_strand");
			break;
		case 's':
			if ((info = calloc(42 + strlen(transcript->ref_id), sizeof(char))) == NULL) {
              			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			sprintf(info, "intronic overlap with %s on opposite strand", transcript->ref_id);
			strcpy(transcript->category, "non-canonical_strand");
			break;
		case 'p':
			if ((info = calloc(42 + strlen(transcript->ref_id), sizeof(char))) == NULL) {
              			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			sprintf(info, "downstream of %s - possible run-on fragment", transcript->ref_id);
			strcpy(transcript->category, "run-on");
			break;
		case 'i':
			if ((info = calloc(25 + strlen(transcript->ref_id), sizeof(char))) == NULL) {
              			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			sprintf(info, "intronic - derived from %s", transcript->ref_id);
			strcpy(transcript->category, "intronic");
			break;
		case 'y':
			if ((info = calloc(31 + strlen(transcript->ref_id), sizeof(char))) == NULL) {
              			printf("memory allocation error in get_metadata()\n");
        			exit(0);
			}
			sprintf(info, "contains %s within its intron(s)", transcript->ref_id);
			strcpy(transcript->category, "wrap-around");
			break;
		case 'c':
		case '=':
		case 'j':
		case 'e':
		case 'o':
		case 'r':
		case 'k':
		case 'm':
		case 'n':
			if (!embl.error) {
				if (embl.start_exon_error) {
					if ((info = calloc(41, sizeof(char))) == NULL) {
              					printf("memory allocation error in get_metadata()\n");
        					exit(0);
					}
					strcpy(info, "ensembl translation start site not found");	
				} else {
					info = calculate_info(transcript, embl);
				}
			} else {
				if ((info = calloc(23, sizeof(char))) == NULL) {
              				printf("memory allocation error in get_metadata()\n");
        				exit(0);
				}
				strcpy(info, "ensembl info not found");
			}
			break;
		default:
			printf("Unknown class_code %c for transcript %s", transcript->class_code, transcript->transcript_id);
			exit(0);
	}
	if (strlen(info) > 511) {
		printf("info too big to store in transcript metadata: %d\n", (int)strlen(info));
		exit(0);
	}
	strcpy(transcript->metadata, info);
	free(info);
	
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
Region find_region(int target, embl_info ref)
{
	int i, flag = 0;

	if (ref.strand == UNKNOWN) 
		return NOT_FOUND;
	if (ref.strand == PLUS) {
		if (target < ref.exons[0].st)
			return UPSTREAM;
		if (target > ref.exons[ref.exon_cnt-1].en)
			return DOWNSTREAM;
	}
	else {
		if (target > ref.exons[0].en)
			return UPSTREAM;
		if (target < ref.exons[ref.exon_cnt-1].st)
			return DOWNSTREAM;
	}
//if not in exon, must be intron
	for (i = 0; i < ref.exon_cnt; ++i)
		if ((target < ref.exons[i].en +1) && (target > ref.exons[i].st -1))
			flag = 1;
	if (!flag)
		return INTRON;
	if (!(ref.translation.st))
		return NONCODING;
	if (ref.strand == PLUS) {
		if (target < ref.translation.st)
			return UTR5;
		if (target > ref.translation.en)
			return UTR3;
	}
	else {
		if (target > ref.translation.en)
			return UTR5;
		if (target < ref.translation.st)
			return UTR3;
	}
	if ((target < ref.translation.en +1) && (target > ref.translation.st -1))
		return TRANSLATED;

	return NOT_FOUND;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int find_frame(int target, embl_info ref)
{
	int i, frame = -1, exon_index = -1, bfs = 0;

	if (ref.strand == PLUS) {
		if (target < ref.exons[0].st)
			exon_index = 0;
		else {
			i = 0;
			while (i < ref.exon_cnt && target > ref.exons[i].st -1) 
				i++;
			exon_index = i - 1;
		}
		bfs = ref.exons[exon_index].bases_from_start + (target - ref.exons[exon_index].st);
	}
	else if (ref.strand == MINUS) {
		if (target > ref.exons[0].en)
			exon_index = 0;
		else {
			i = 0;
			while (i < ref.exon_cnt && target < ref.exons[i].en +1) 
				i++;
			exon_index = i - 1;
		}		
		bfs = ref.exons[exon_index].bases_from_start + (ref.exons[exon_index].en - target);
	}
	else {
		printf("Problem in find_frame() - strand unknown\n");
		exit(0);
	}

	if (bfs > -1)
		frame = bfs%3;
	else {
		bfs = bfs%3;
		frame = (bfs == 0) ? 0 : 3 + bfs;
	}
	if (frame < 0 || frame > 2) {
		printf("Frame is not right - %d\nbfs %d. Transcript %s.\n", frame, bfs, ref.id);
		exit(0);
	}
	return frame;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
char *calculate_info(genome_info *peptide, embl_info ref)
{
	char tmp[512], *info = NULL, chunk[512], full_id[50];
	Region start, end;
	int frame = -1;

	full_id[0] = '\0';
	strcpy(full_id, ref.id);
	if (ref.version) {
		tmp[0] = '\0';
		sprintf(tmp, ".%d", ref.version);
		strcat(full_id, tmp);
	}
// 	report in info string - peptide found in %s region or upstream/downstream of ensembl transcript, in-frame/frame 1/frame 2 with reference to ensembl transcript (where %s is translated, 5'UTR, 3'UTR)  

	tmp[0] = chunk[0] = '\0';
//check strand == strand
	if (peptide->strand != ref.strand) {
		printf("cuffcompare reference transcript strand does not match the Ensembl transcript strand: possible incompatibility with latest Ensembl version\n");
		printf("cuffcompare tracking file %s, Ensembl reference %s\n", peptide->transcript_id, ref.id);
		printf("pep strand %d, embl strand %d\n", peptide->strand, ref.strand);
		strcpy(peptide->category, "error");
		strcpy(tmp, "Cuffcompare reference transcript does not match Ensembl transcript info. Consider manual investigation");
		if ((info = calloc((strlen(tmp) + 1), sizeof(char))) == NULL) {
        		printf("memory allocation error in calculate_info()\n");
        		exit(0);
		}
		strcpy(info, tmp);
		return info;
	}
	
//find whether peptide start and end is in translation chunk, utr chunk or outside ref transcript etc
	start = find_region(peptide->start, ref);
	end = find_region(peptide->end, ref);

	if (start == end) {
		switch(start) {
			case UPSTREAM:
				sprintf(tmp, "peptide found upstream of %s", full_id); 
				strcpy(peptide->category, "upstream");
				break;
			case DOWNSTREAM:
				sprintf(tmp, "peptide found downstream of %s", full_id); 
				strcpy(peptide->category, "downstream");
				break;
			case UTR5:
				sprintf(tmp, "peptide found in 5'-UTR of %s", full_id); 
				strcpy(peptide->category, "5'UTR");
				break;
			case UTR3:
				sprintf(tmp, "peptide found in 3'UTR of %s", full_id); 
				strcpy(peptide->category, "3'UTR");
				break;
			case TRANSLATED:
				sprintf(tmp, "peptide found in translated region of %s", full_id); 
				strcpy(peptide->category, "translated");
				break;
			case INTRON:
				sprintf(tmp, "peptide found in intronic region of %s", full_id); 
				strcpy(peptide->category, "intron");
				break;
			case NONCODING:
				sprintf(tmp, "peptide found in non-coding exon of %s", full_id); 
				strcpy(peptide->category, "non-coding");
				break;
			case NOT_FOUND:
				printf("Something's wrong with finding region\n");
				exit(0);
		}	
	}
	else {
		strcpy(peptide->category, "junction");
		strcpy(tmp, "Peptide starts ");
		switch(start) {
			case UPSTREAM:
				strcat(tmp, "upstream and ends "); 
				break;
			case DOWNSTREAM:
				strcat(tmp, "downstream and ends "); 
				break;
			case UTR5:
				strcat(tmp, "in 5'-UTR and ends "); 
				break;
			case UTR3:
				strcat(tmp, "in 3'-UTR and ends "); 
				break;
			case TRANSLATED:
				strcat(tmp, "in translated region and ends "); 
				break;
			case INTRON:
				strcat(tmp, "in intronic region and ends "); 
				break;
			case NONCODING:
				strcat(tmp, "in non-coding exon and ends "); 
				break;
			case NOT_FOUND:
				printf("Something's wrong with finding region\n");
				exit(0);
		}	
		switch(end) {
			case UPSTREAM:
				sprintf(chunk, "upstream of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case DOWNSTREAM:
				sprintf(chunk, "downstream of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case UTR5:
				sprintf(chunk, "in 5'-UTR of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case UTR3:
				sprintf(chunk, "in 3'-UTR of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case TRANSLATED:
				sprintf(chunk, "in translated region of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case INTRON:
				sprintf(chunk, "in intronic region of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case NONCODING:
				sprintf(chunk, "in non-coding exon of %s", full_id);
				strcat(tmp, chunk); 
				break;
			case NOT_FOUND:
				printf("Something's wrong with finding region\n");
				exit(0);
		}	
		
	}
	
//find whether peptide start is in frame with ref transcript
	if (ref.translation.embl_st) {
		if ((frame = find_frame(peptide->start, ref)) < 0) {
			printf("Can't find frame...\n");
			exit(0);
		}
		if (frame >2) {
			printf("frame %d, pepref %s, transcript %s error? %d\n", frame, peptide->ref_id, ref.id, ref.error);
		}
		switch(frame) {
			case(0):
				strcat(tmp, ", in-frame with canonical translation");
				strcat(peptide->category, "_in-frame");
				break;
			case(1):
				strcat(tmp, ", frame +1 relative to canonical translation"); 
				strcat(peptide->category, "_out-of-frame");
				break;
			case(2):
				strcat(tmp, ", frame +2 relative to canonical translation"); 
				strcat(peptide->category, "_out-of-frame");
				break;
			default:
				printf("Problem with find_frame(); frame = %d\n", frame);
				exit(0);
		}	
	}

	if ((info = calloc((strlen(tmp) + 1), sizeof(char))) == NULL) {
        	printf("memory allocation error in calculate_info()\n");
        	exit(0);
	}
	strcpy(info, tmp);
	return info;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int get_transcripts(FILE *f, seq_info **seqs, int cnt, int file_num, transcriptome_mutations *mut_from_vcf)
{
	char **aa_seq = NULL, transcript[MAXSEQLEN], line[MAXSTR], title[MAXSTR], *pnt = NULL, *ptr = NULL, copy[MAXSEQLEN/3];
	int i, j, len = 0, first = 1, tlen = 0, seq_cnt = 0, pos = 0;	
	genome_info *tmp;

        transcript[0] = title[0] = copy[0] = '\0';
	
	if ((aa_seq = calloc(3, sizeof(char *))) == NULL) {
              	printf("memory allocation error in get_transcripts()\n");
        	exit(0);
	}
	
        while (fgets(line, MAXSTR-1, f) != NULL) {
                len = (int)strlen(line);
                while (line[len - 1] == '\n' || line[len - 1] == '\r') {
                        line[len -1] = '\0';
                        len = (int)strlen(line);
                }
		
                if (line[0] == '>') {
                        if (first) {
                                strcpy(title, line +1);
                                first = 0;
                        }
			else {
				seq_cnt ++;
				for (i = 0; i < 3; i++) {
                                	if ((aa_seq[i] = calloc(((tlen+1)/3 + 1), sizeof(char))) == NULL) {
                                        	printf("memory allocation error in get_transcripts()\n");
                                                exit(0);
                                        }
				}
	//make translation, check for match, store info including coordinates!, then save new title and reset
				make_translation(&aa_seq, transcript);
				for (i = 0; i < cnt; ++i) { //for each peptide
					for (j = 0; j < 3; j++) { //for each frame translation
						copy[0] = '\0'; 
						strcpy(copy, aa_seq[j]);
						ptr = copy;
                                		while ((pnt = (strstr(ptr, (*seqs)[i].no_ptm))) != NULL) { //look for peptide in translation
							pos = ((strlen(aa_seq[j]) - strlen(pnt)) * 3) + j;
							(*seqs)[i].transcriptomes[file_num].transcript_cnt++;
                                                	if ((*seqs)[i].transcriptomes[file_num].transcript_cnt == 1) {
                                                        	if (((*seqs)[i].transcriptomes[file_num].transcripts = calloc(1, sizeof(genome_info))) == NULL) {
                                                        	        printf("memory allocation error in get_transcripts()\n");
                                                        	        exit(0);
                                                        	}
                                                	}
                                                	else {
                                                        	if ((tmp = realloc((*seqs)[i].transcriptomes[file_num].transcripts, (((*seqs)[i].transcriptomes[file_num].transcript_cnt) * sizeof(genome_info)))) == NULL) {
                                                                	printf("memory allocation error in get_accessions()\n");
                                                                	exit(0);
                                                        	}
                                                        	(*seqs)[i].transcriptomes[file_num].transcripts = tmp;
								memset((*seqs)[i].transcriptomes[file_num].transcripts + (*seqs)[i].transcriptomes[file_num].transcript_cnt - 1, 0, sizeof(genome_info));	
                                                	}
							parse_title(&((*seqs)[i].transcriptomes[file_num].transcripts[(*seqs)[i].transcriptomes[file_num].transcript_cnt - 1]), title, pos, pos -1 + (strlen((*seqs)[i].no_ptm)*3), mut_from_vcf);
							ptr = pnt + 1;
						}
					}
				}
                        	title[0] = '\0';
				for (i = 0; i < 3; i++) {
					if (aa_seq[i] != NULL)
						free(aa_seq[i]);
				}
                        	strcpy(title, line +1);
                        	tlen=0;
                        	transcript[0]='\0';
			}
		}
                else  {
                        tlen += len;
                        if (tlen > MAXSEQLEN) {
                                printf("Increase MAXSEQLEN\n");
                                exit(0);
                        }
                        strcat(transcript, line);
                }
	}
//deal with last!	
	//make translation, check for match, store,
	seq_cnt ++;
	for (i = 0; i < 3; i++) {
               	if ((aa_seq[i] = calloc(((tlen+1)/3 + 1), sizeof(char))) == NULL) {
                       	printf("memory allocation error in get_transcripts()\n");
                        exit(0);
                }
	}
	make_translation(&aa_seq, transcript);
	for (i = 0; i < cnt; ++i) { //for each peptide
		for (j = 0; j < 3; j++) { //for each frame translation
			copy[0] = '\0'; 
			strcpy(copy, aa_seq[j]);
                       		while ((pnt = (strstr(copy, (*seqs)[i].no_ptm))) != NULL) { //look for peptide in translation
					pos = ((strlen(aa_seq[j]) - strlen(pnt)) * 3) + j;
					(*seqs)[i].transcriptomes[file_num].transcript_cnt++;
                                        if ((*seqs)[i].transcriptomes[file_num].transcript_cnt == 1) {
                                        	if (((*seqs)[i].transcriptomes[file_num].transcripts = calloc(1, sizeof(genome_info))) == NULL) {
                                                	printf("memory allocation error in get_transcripts()\n");
                                                        exit(0);
                                                }
                                        }
                                        else {
                                              	if ((tmp = realloc((*seqs)[i].transcriptomes[file_num].transcripts, (((*seqs)[i].transcriptomes[file_num].transcript_cnt) * sizeof(genome_info)))) == NULL) {
                                                       	printf("memory allocation error in get_accessions()\n");
                                                       	exit(0);
                                               	}
                                               	(*seqs)[i].transcriptomes[file_num].transcripts = tmp;
						memset((*seqs)[i].transcriptomes[file_num].transcripts + (*seqs)[i].transcriptomes[file_num].transcript_cnt - 1, 0, sizeof(genome_info));	
                                       	}
					parse_title(&((*seqs)[i].transcriptomes[file_num].transcripts[(*seqs)[i].transcriptomes[file_num].transcript_cnt - 1]), title, pos, pos -1 + (strlen((*seqs)[i].no_ptm)*3),mut_from_vcf);
					copy[0] = '\0';
					strcpy(copy, pnt + 1); 
				}
		}
	}
	
	for (i = 0; i < 3; i++) {
		if (aa_seq[i] != NULL)
			free(aa_seq[i]);
	}
	if (aa_seq != NULL) {
		free(aa_seq);
	}
	return seq_cnt;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void parse_title(genome_info *transcript, char *title, int st, int en, transcriptome_mutations *mut_from_vcf)
{
	exon *exons = NULL;
        char *ptr = NULL, *pnt = NULL, chunk[strlen(title)], chunk2[strlen(title)], c = '\0', loc[50];
        int i, j, exon_cnt = 0, flag = 0, filtered_cnt=0, match = -1;
	vcf_mut_info *filtered_mut;

	chunk[0] = chunk2[0] = '\0';	
//get id and loc
	sscanf(title, "%s %s", transcript->transcript_id, loc);
//get strand
        c = loc[strlen(loc)-1];
        if (c == '+')
                transcript->strand = PLUS;
        else if (c == '-')
                transcript->strand = MINUS;
        else {
                transcript->strand = UNKNOWN;
        }
//get chr
        ptr = loc + 4;
        pnt = move_ptr(ptr, '|');
        strcpy(transcript->chr, ptr);
        *pnt = '|';
        chunk[0] = chunk2[0] = '\0';
//get exons
	if ((ptr = strstr(title, "exons:")) != NULL) {
        	ptr += 6;
                pnt = move_ptr(ptr, ' ');
                strcpy(chunk, ptr);
                *pnt = ' ';
                ptr = pnt + 6;
		if (strchr(ptr,' ') != NULL )
                	pnt = move_ptr(ptr,' ');
                strcpy(chunk2, ptr);
                *pnt = ' ';
                exon_cnt = parse_exons(&exons, chunk, chunk2, transcript->strand);
		if (exon_cnt < 1) {
                	printf("parse_exons failed...\n");
                        exit(0);
                }
//check mutations are in exons if mut_from_vcf != NULL
		if (mut_from_vcf != NULL) {
	//find matching contig
			for (i = 0; i < mut_from_vcf->contig_cnt; ++i) {
                		if (strcmp(mut_from_vcf->contigs[i].id, transcript->chr) == 0) {
                        		match = i;
                        	}
               		}
               	 	if (match < 0) {
                		printf("chromosome names don't match vcf in parse_title: chr %s\n", transcript->chr);
                                exit(0);
			}
                	filtered_mut = filter_mut_by_exons(exons, exon_cnt, &(mut_from_vcf->contigs[match]), &filtered_cnt);
		}
	}
//convert transcript coordinates to genome coordinates
	transcript->start = transcript_to_genome_loc(exons, exon_cnt, transcript->strand, st + 1);
	transcript->end = transcript_to_genome_loc(exons, exon_cnt, transcript->strand, en + 1);
	
	if (transcript->start < 0 || transcript->end < 0) {
		printf("Problem converting transcript to genome coordinates: start %d, end %d\n", transcript->start, transcript->end);
		exit(0);
	}
//check if mutation within peptide coordinates and if so copy to 'transcript'
	if (filtered_cnt) {
		if (filtered_mut == NULL) {
			printf("Sanity check: filter_mut_by_exons() didn't work\n");
			exit(0);
		}
		add_mutation_info(transcript, &filtered_mut, filtered_cnt);
	}
	if (exons != NULL)
                free(exons);
	if (filtered_cnt) {
		if (filtered_mut != NULL)
			free(filtered_mut);
	}
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int transcript_to_genome_loc(exon *exons, int exon_cnt, int strand, int pos)
{
	int i, loc = -1;

	for (i = 0; i < exon_cnt; i++) {
		if ((pos < (exons[i].transcript_en +1)) && (pos > (exons[i].transcript_st - 1))) {
			if (strand == PLUS || strand == UNKNOWN) {
				loc = exons[i].g_en - exons[i].transcript_en + pos;
			} else if (strand == MINUS) { 
				loc = exons[i].transcript_en + exons[i].g_st - pos;
                       	} 
		}
	}
	return loc;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
vcf_mut_info *filter_mut_by_exons(exon *exons, int exon_cnt, contig_info *muts_from_vcf, int *filtered_cnt)
{
        int i, j;
	vcf_mut_info *filtered_mut = NULL;

	*filtered_cnt = 0;
//filtered_cnt should work below, because exons shouldn't overlap, so mut should fall only in one exon each
        for (i=0; i<muts_from_vcf->mut_cnt; ++i) {
                for (j=0; j<exon_cnt; ++j) { 
                        if ((muts_from_vcf->mut[i].alt_pos < (exons[j].g_en +1)) && (muts_from_vcf->mut[i].alt_pos > (exons[j].g_st -1))) { 
                        	muts_from_vcf->mut[i].flag = 1;
				*filtered_cnt += 1;
			}
                }
        }
	if (*filtered_cnt) {
        	if ((filtered_mut=calloc(*filtered_cnt, sizeof(vcf_mut_info)))==NULL) {
                	printf("memory allocation error in filter_mut_by_exons\n");
                	exit(0);
        	}
		j=0;
        	for (i=0; i<muts_from_vcf->mut_cnt; ++i) {
			if (muts_from_vcf->mut[i].flag) {
				if (j > *filtered_cnt + 1) {
					printf("Greater than one exon per mutation?\n");
					exit(0);
				}	
				filtered_mut[j].alt_pos = muts_from_vcf->mut[i].alt_pos;
				strcpy(filtered_mut[j].mutation, muts_from_vcf->mut[i].mutation);
				j++;
			}
		}
		reset_flags(muts_from_vcf);
	}
        return filtered_mut;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
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
                (*exons)[i].g_st = atoi(ptr);
                ptr = pnt +1;
                pnt = move_ptr(ptr, ',');
                (*exons)[i].g_en = atoi(ptr);
                if (i < cnt -1)
                        ptr = pnt +1;
        }
        ptr = chunk2;
//copy transcript loci into exon structure array
        if (strand == PLUS || strand == UNKNOWN) {
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
                (*exons)[i].len = (*exons)[i].g_en - (*exons)[i].g_st + 1;
                if ((*exons)[i].len != ((*exons)[i].transcript_en - (*exons)[i].transcript_st + 1)) {
                        printf("Exon and seg length don't match!\n");
                        exit(0);
        	}
        }
        return cnt;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void make_translation(char ***aa, char *transcript)
{
	int x,y,len, tlen = strlen(transcript);
	char c[5][5][5] = {{"FFLLX","SSSSX","YYZZX","CCZWX","XXXXX"},
                   {"LLLLX","PPPPX","HHQQX","RRRRX","XXXXX"},
                   {"IIIMX","TTTTX","NNKKX","SSRRX","XXXXX"},
                   {"VVVVX","AAAAX","DDEEX","GGGGX","XXXXX"},
                   {"XXXXX","XXXXX","XXXXX","XXXXX","XXXXX"}};
	int coded[] = {2,4,1,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,0,0,4,4,4,4,4};

	for(x = 0; x < 3; x++){
        	len = 0;
        	for(y = x;y < tlen; y += 3) {
                	if (transcript[y] != '\0' && transcript[y + 1] != '\0' && transcript[y + 2] != '\0') {
	                        (*aa)[x][len] = c[coded[transcript[y] - 'A']][coded[transcript[y + 1] - 'A']][coded[transcript[y + 2] - 'A']];
		                len ++;
                        	(*aa)[x][len] = '\0';
                        }
                }
        }
}
/**********************************************************************************************************/

/**********************************************************************************************************/
void get_uniprot(FILE *f, seq_info **seqs, int cnt)
{
        int i, len = 0, first = 1, protlen = 0;
        char line[MAXSTR], protein[MAXPROT], title[512], **tmp = NULL;

        protein[0] = title[0] = '\0';

        while (fgets(line, MAXSTR-1, f) != NULL) {
                len = (int)strlen(line);
                while (line[len - 1] == '\n' || line[len - 1] == '\r') {
                        line[len -1] = '\0';
                        len = (int)strlen(line);
                }
                if (line[0] == '>') {
                        if (first) {
                                strcpy(title, line +1);
                                first = 0;
                        }
                        else {
                                for (i = 0; i < cnt; i++) {
                                	if ((strstr(protein, (*seqs)[i].no_ptm)) != NULL) {
                                        	(*seqs)[i].uniprot_cnt++;
                                                if ((*seqs)[i].uniprot_cnt == 1) {
                                                	if (((*seqs)[i].uniprot_accessions = calloc(1, sizeof(char *))) == NULL) {
                                                        	printf("memory allocation error in get_accessions()\n");
                                                                exit(0);
                                                        }
                                               	}
                                                else {
                                                	if ((tmp = realloc((*seqs)[i].uniprot_accessions, (((*seqs)[i].uniprot_cnt) * sizeof(char *)))) == NULL) {
                                                        	printf("memory allocation error in get_accessions()\n");
                                                                exit(0);
                                                        }
                                                        (*seqs)[i].uniprot_accessions = tmp;
						}
                                                if (((*seqs)[i].uniprot_accessions[(*seqs)[i].uniprot_cnt-1] = calloc((int)(strlen(title) + 1), sizeof(char))) == NULL) {
                                                	printf("memory allocation error in get_accessions()\n");
                                                        exit(0);
                                                }
                                                strcpy((*seqs)[i].uniprot_accessions[(*seqs)[i].uniprot_cnt -1], title);
                                       	}
                             	}
                                title[0] = '\0';
                                strcpy(title, line +1);
                                protlen=0;
                                protein[0]='\0';
                        }
                }
                else  {
                        protlen += len;
                        if (protlen > MAXPROT) {
                                printf("Increase MAXPROT\n");
                                exit(0);
                        }
                        strcat(protein, line);
                }
        }
//last protein
	for (i = 0; i < cnt; i++) {
        	if ((strstr(protein, (*seqs)[i].no_ptm)) != NULL) {
                	(*seqs)[i].uniprot_cnt++;
                        if ((*seqs)[i].uniprot_cnt == 1) {
                        	if (((*seqs)[i].uniprot_accessions = calloc(1, sizeof(char *))) == NULL) {
                                	printf("memory allocation error in get_accession()\n");
                                        exit(0);
                                }
                        }
                        else {
                        	if ((tmp = realloc((*seqs)[i].uniprot_accessions, (((*seqs)[i].uniprot_cnt) * sizeof(char *)))) == NULL) {
                                	printf("memory allocation error in get_accession()\n");
                                        exit(0);
                                }
                                (*seqs)[i].uniprot_accessions = tmp;
                       	}
                        if (((*seqs)[i].uniprot_accessions[(*seqs)[i].uniprot_cnt-1] = calloc((int)(strlen(title) + 1), sizeof(char))) == NULL) {
                        	printf("memory allocation error in get_accession()\n");
                                exit(0);
                        }
                        strcpy((*seqs)[i].uniprot_accessions[(*seqs)[i].uniprot_cnt -1], title);
                }
        }
	return;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int sort_alpha(const void *a, const void *b) {
        seq_info *e1 = (seq_info *)a, *e2 = (seq_info *)b;
    return strcmp(e1->seq, e2->seq);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
char *move_ptr(char *str,char delim)
{
char *tmp = NULL;

if((tmp = strchr(str,delim)) != NULL) *tmp = '\0';

return(tmp);
}
/**********************************************************************************************************/

/**********************************************************************************************************/
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
/**********************************************************************************************************/

/**********************************************************************************************************/
vcf_info *read_alt_vcf(FILE *f, contig_info **contigs, int lines_cnt, int contig_cnt)
{
        vcf_info *tmp = NULL;
        int i = 0, j = 0, k=0;
        char line[MAXSTR], *pnt = NULL, *ptr = NULL;

        if ((*contigs=(contig_info *)calloc(contig_cnt, sizeof(contig_info)))==NULL) {
                printf("memory allocation error in read_alt_vcf\n");
                exit(0);
        }
        for (k = 0; k<contig_cnt; ++k)
                (*contigs)[k].mut_cnt = 0;

        if ((tmp=(vcf_info *)calloc(lines_cnt, sizeof(vcf_info)))==NULL) {
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
/**********************************************************************************************************/

/**********************************************************************************************************/
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
			contigs[i].mut[j].alt_pos = v[k].pos + offset;
                        sprintf(contigs[i].mut[j].mutation, "%s->%s@%d(ref:%d)", v[k].ref, v[k].alt, (v[k].pos + offset), v[k].pos);
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
/**********************************************************************************************************/

/**********************************************************************************************************/
void add_deleted_pos_memory(int **del, int old, int new)
{
        int *tmp = NULL;

        if (*del == NULL) {
                if ((*del = (int *)calloc(new, sizeof(int))) == NULL) {
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
/**********************************************************************************************************/

/**********************************************************************************************************/
static size_t receive(void *contents, size_t size, size_t nmemb, void *data)
{
        size_t realsize = size * nmemb;

        curl_data *mem = (curl_data *)data;

        char *ptr = realloc(mem->memory, mem->size + realsize + 1);
        if(ptr == NULL) {
    /* out of memory! */
                printf("not enough memory (realloc returned NULL)\n");
                return 0;
        }

        mem->memory = ptr;
        memcpy(&(mem->memory[mem->size]), contents, realsize);
        mem->size += realsize;
        mem->memory[mem->size] = 0;

        return realsize;
}
/**********************************************************************************************************/

/**********************************************************************************************************/
int liftover(int old, contig_info contig)
{
	int new = 0, index = 0;
	
	while ((contig.mut[index].pos < old) && index < contig.mut_cnt) {
        	index++;
        }
	if (index != 0) {
        	new = old + contig.mut[index-1].offset;
        }

	return new;
}
/**********************************************************************************************************/
