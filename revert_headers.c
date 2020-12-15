#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#define MAXSTR 51200

int count_headers(FILE *);

/* revert_headers: FastaAlternateReferenceMaker changes the chromosome names; this changes
them back to match the reference genome */
/*********************************************************************************************/
void main(int argc,char **argv)
{
        FILE *f, *g;
	char **header = NULL, line[MAXSTR];
	int header_cnt = 0, cnt = 0, i = 0;

	if (argc < 3) {
		printf("Specify reference genome then new genome you want to fix\n");
		exit(0);
	}
	printf("Reverting headers from %s to headers from %s\n", argv[2], argv[1]);
	if ((f = fopen(argv[1],"r")) == NULL) {
                printf("Can't open file %s\n",argv[1]);
                exit(0);
        }
	if ((g = fopen("tmpc.fasta","w")) == NULL) {
                printf("Can't open file tmpc.fasta\n");
                exit(0);
        }
	
	header_cnt = count_headers(f);
	printf("header_cnt = %d\n", header_cnt);
	if ((header=calloc(header_cnt, sizeof(char *)))==NULL) {
                printf("memory allocation error in main\n");
                exit(0);
        }
	i = 0;
	while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>' && i < header_cnt) {
			if ((header[i]=calloc((strlen(line) + 1), sizeof(char)))==NULL) {
                                printf("memory allocation error in main\n");
                                exit(0);
                        }
                        strcpy(header[i], line);
                        i++;
		}
	}
	fclose(f);

	if ((f = fopen(argv[2],"r")) == NULL) {
                printf("Can't open file %s\n",argv[2]);
                exit(0);
        }
	i=0;
	while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>') {
			cnt++;
			if (cnt > header_cnt) {
				printf("the number of headers in the two genomes don't match\n");
				exit(0);
			}
			fprintf(g, "%s", header[i]);
			printf("Was %snow %s\n", line, header[i]);
			i++;
		}
		else
			fprintf(g, "%s",line);
	}
	printf("Finished!\n");
	for (i = 0; i < header_cnt; i++) {
		if (header[i] != NULL) 
			free(header[i]);
	}
	if (header != NULL) 
		free(header);

	fclose(f);
	fclose(g);
	exit(0);
}
/*********************************************************************************************/

/*********************************************************************************************/
int count_headers(FILE *f)
{
        int cnt = 0;
        char line[MAXSTR];

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>')
                        cnt+=1;
        }
        rewind(f);

        return(cnt);
}
/*********************************************************************************************/
