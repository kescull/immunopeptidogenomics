#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#define MAXSTR 10000
#define TOLERANCE 0.1
#define MAXIONS 2500
#define MININT 0.01
#define PEAK 3
#define MAYBE 1
#define MATCH 2
#define NONE 0
#define COMMON 1
#define ALL 0

typedef struct ions ions;
struct ions {
	double mz;
	double intensity;
	int tag;
	}; 

ions *find_ions(FILE *, int *);
int make_vectors(FILE *,ions *, int, ions *, int, double **, double **, int *,int);
int sort_by_mz(const void *,const void *);
void check_peak(ions *, int); 
void check_match(ions *, int, ions *, int, int *);
void make_match(ions *, ions *, ions *, ions *, int *);
double dotproduct(double *, double *, int);
ions *clean_spectrum(ions *, int, int *);

/* msDot.c: compare pairs of dta file MS/MS spectra using dot product. 
Input: Directory containing pairs dta files with naming convention *_mhc_peptide.dta and *_synth_peptide.dta,
where * is an identifier such as the presumed peptide sequence. 
Program reads ions, cleans spectra to remove peaks with intensity below MININT * max intensity,
and makes vectors with equal numbers of values in each array (ie. same number of peaks), 
matching peaks within TOLERANCE between spectra. Calculates normalised dot product of the vectors, ensuring
Normalised dot product = a∙b / (|a||b|) where |a| = √(a∙a) and |b| = √(b∙b)
Outputs value between 0 and 1 where closer to 1 means more similar.*/

/**********************************************************************************************************/
void main(int argc,char **argv)
{
	DIR *dirp;
        struct dirent *d;
	FILE *f, *g, *h;
	ions *spectrum1 = NULL, *spectrum2 = NULL, *clean_spec1 = NULL, *clean_spec2 = NULL;
	int ion_cnt1=0, ion_cnt2=0, i, match_cnt=0, least_ions = 0, clean_cnt1 = 0, clean_cnt2 = 0, common = 0;
 	double *v1=NULL, *v2=NULL, prop = 0.0;
	double len1 = 0.0, len2 = 0.0, dot = 0.0;
	char peptide[100], fn[512], *ptr, copy[512], output[512];

	if (!(argc==2)) {
                printf("Error: specify directory\n");
                exit(0);
        }

	if ((dirp = opendir(argv[1]))==NULL) {
                printf("error: check directory\n");
                exit(0);
        }
	sprintf(output,"%s/dotproducts_full.csv", argv[1]);
	if((g = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n",output);
                exit(0);
        }
	sprintf(output,"%s/dotproducts_only.csv", argv[1]);
	if((h = fopen(output,"w")) == NULL) {
                printf("Can't open file %s\n",output);
                exit(0);
        }
	fprintf(h,"Dotproducts of normalised spectra\nPeptide,Dotproduct,Common/least ions\n");

        while ((d = readdir(dirp)) != NULL) {
                if((strstr(d->d_name,"_mhc_peptide.dta")) != NULL) {
                        //get spectrum1
			ion_cnt1 = ion_cnt2 = least_ions = clean_cnt1 = clean_cnt2 = common = 0;
			sprintf(fn, "%s/%s", argv[1], d->d_name);
			if((f = fopen(fn,"r")) == NULL) {
                		printf("Can't open file %s\n",fn);
                		exit(0);
			}
			spectrum1 = find_ions(f, &ion_cnt1);

			fclose(f);
			
			//get spectrum2
			copy[0] = '\0';
			fn[0] = '\0';
			strcpy(copy, d->d_name);
			ptr = strstr(copy,"_mhc_peptide.dta");
			*ptr = '\0';
			strcpy(peptide,copy);
			sprintf(fn, "%s/%s%s", argv[1], peptide, "_synth_peptide.dta");

			if((f = fopen(fn,"r")) == NULL) {
                		printf("Can't open file %s\n",fn);
                		exit(0);
			}
			spectrum2 = find_ions(f, &ion_cnt2);
	
			fclose(f);

			fprintf(g,"%s\n", peptide);
			fprintf(g,"MHC spectrum ion count = %d\n", ion_cnt1);
			fprintf(g,"Synthetic spectrum ion count = %d\n", ion_cnt2);
			
			clean_spec1 = clean_spectrum(spectrum1,ion_cnt1,&clean_cnt1);
			clean_spec2 = clean_spectrum(spectrum2,ion_cnt2,&clean_cnt2);

			least_ions = clean_cnt1 < clean_cnt2 ? clean_cnt1 : clean_cnt2;
			fprintf(g,"After:\nclean MHC spectrum ion count = %d\n", clean_cnt1);
			fprintf(g,"clean Synthetic spectrum ion count = %d\n", clean_cnt2);
			fprintf(g,"Smallest spectrum = %d ions\n", least_ions);

			qsort(clean_spec1, clean_cnt1, sizeof(ions), sort_by_mz);
			qsort(clean_spec2, clean_cnt2, sizeof(ions), sort_by_mz);

			match_cnt = make_vectors(g, clean_spec1, clean_cnt1, clean_spec2, clean_cnt2, &v1, &v2,&common,0);
			fprintf(g,"Common %d\n", common);
			prop = (double)common/least_ions;
			len1 = sqrt(dotproduct(v1, v1, match_cnt));
			len2 = sqrt(dotproduct(v2, v2, match_cnt));
			dot = dotproduct(v1, v2, match_cnt)/(len1 * len2);
			fprintf(g,"Dotproduct of normalised spectra for %s = %lf\n", peptide, dot);
			fprintf(h,"%s,%lf,%lf\n", peptide, dot, prop);
			
			//clean up
			ion_cnt1 = ion_cnt1 = 0;
			if (spectrum1 != NULL)
				free(spectrum1);
			if (spectrum2 != NULL)
				free(spectrum2);
			if (clean_spec1 != NULL)
				free(clean_spec1);
			if (clean_spec2 != NULL)
				free(clean_spec2);
			if(v1 != NULL)
				free(v1);	
			if(v2 != NULL)
				free(v2);	
		}
	}
	fclose(g);
	fclose(h);
        closedir(dirp);
	exit(0);
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
ions *clean_spectrum(ions *spectrum, int cnt, int *clean_cnt)
{
        double max = 0.0;
        int i;
        ions tmp[cnt];
        ions *new = NULL;

        *clean_cnt=0;
        for (i=0; i<cnt; ++i) {
                if (spectrum[i].intensity > max)
                        max = spectrum[i].intensity;
        }
	printf("Max intensity: %lf\n", max);
        for (i=0; i<cnt; ++i) {
                if (spectrum[i].intensity > (MININT * max)) {
                        tmp[*clean_cnt] = spectrum[i];
                        *clean_cnt +=1;
                }
                else {
                        printf("Fail Ion %d) %lf %lf\n", i, spectrum[i].mz, spectrum[i].intensity);
                }
        }
        if (*clean_cnt>0) {
                if((new = (ions *)calloc(*clean_cnt,sizeof(ions))) == NULL) {
                        printf("Memory allocation error\n");
                        exit(0);
}
                memcpy(new, tmp, (*clean_cnt * sizeof(ions)));
        }
        return new;
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
ions *find_ions(FILE *f, int *ion_cnt)
{
	ions tmp[MAXIONS], *spectrum = NULL;
	char line[MAXSTR];
	int first = 1;

	while(fgets(line,MAXSTR - 1, f) != NULL) {
		if (first) {
			first = 0;
		} else {
			sscanf(line,"%lf %lf",&tmp[*ion_cnt].mz, &tmp[*ion_cnt].intensity);
			tmp[*ion_cnt].tag = NONE;
                        *ion_cnt+=1;

                        if(*ion_cnt >= MAXIONS) {
                        	printf("Error - increase the size of MAXIONS (it's %d)\n", *ion_cnt);
                        	exit(0);
                        }
             	}
	}
	if((spectrum = calloc(*ion_cnt,sizeof(ions))) == NULL) {
        	printf("Memory allocation error\n");
        	exit(0);
        }
        memcpy(spectrum, tmp, sizeof(ions) * (*ion_cnt));
	return (spectrum);
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
int make_vectors(FILE *g, ions *spectrum1, int ion_cnt1, ions *spectrum2, int ion_cnt2, double **v1, double **v2, int *common, int mode)
{
	int i=0, j=0, k=0, l=0, match_cnt=0, remaining=0;
	ions matches1[(ion_cnt1+ion_cnt2)], matches2[(ion_cnt1+ion_cnt2)];

	check_match(spectrum1, ion_cnt1, spectrum2, ion_cnt2, &remaining); /*initiate tag as 'none' then label possible matches 'maybe'*/
	check_peak(spectrum1, ion_cnt1); /*check 'maybe' ion intensities for those at least twice as big as their neighbours*/ 
	check_peak(spectrum2, ion_cnt2);

	for (i=0; i<ion_cnt1; i++) { /*match peaks to peaks - and surrounding*/
                for (j=0; j<ion_cnt2; j++) {
			if ((spectrum1[i].tag==PEAK && spectrum2[j].tag==PEAK) && ((isless((spectrum1[i].mz - spectrum2[j].mz), TOLERANCE)) && (isgreater((spectrum1[i].mz - spectrum2[j].mz), -TOLERANCE)))) {
				make_match(&spectrum1[i], &spectrum2[j], matches1, matches2, &match_cnt);
				for (k=i-1, l=j-1; k>-1 && l>-1; --k, --l) {
					if (spectrum1[k].tag == MAYBE && spectrum2[l].tag == MAYBE && (isless((spectrum1[k].mz - spectrum2[l].mz), TOLERANCE) && isgreater((spectrum1[k].mz - spectrum2[l].mz), -TOLERANCE))) 
						make_match(&spectrum1[k], &spectrum2[l], matches1, matches2, &match_cnt);
					else break;
				}
				for (k=i+1, l=j+1; k<(ion_cnt1) && l<(ion_cnt2); ++k, ++l) {
					if (spectrum1[k].tag == MAYBE && spectrum2[l].tag == MAYBE && (isless((spectrum1[k].mz - spectrum2[l].mz), TOLERANCE) && isgreater((spectrum1[k].mz - spectrum2[l].mz), -TOLERANCE)))
						make_match(&spectrum1[k], &spectrum2[l], matches1, matches2, &match_cnt);
					else break;
				}
			}
		}
	}
	for (i=0; i<ion_cnt1; i++) { /*match remaining peaks - and surrounding*/
                for (j=0; j<ion_cnt2; j++) {
			if (((spectrum1[i].tag == PEAK && spectrum2[j].tag == MAYBE) || (spectrum1[i].tag == MAYBE && spectrum2[j].tag == PEAK)) && ((isless((spectrum1[i].mz - spectrum2[j].mz), TOLERANCE)) && (isgreater((spectrum1[i].mz - spectrum2[j].mz), -TOLERANCE)))) {
				make_match(&spectrum1[i], &spectrum2[j], matches1, matches2, &match_cnt);
				for (k=i-1, l=j-1; k>-1 && l>-1; --k, --l) {
					if (spectrum1[k].tag == MAYBE && spectrum2[l].tag == MAYBE && (isless((spectrum1[k].mz - spectrum2[l].mz), TOLERANCE) && isgreater((spectrum1[k].mz - spectrum2[l].mz), -TOLERANCE))) 
						make_match(&spectrum1[k], &spectrum2[l], matches1, matches2, &match_cnt);
					else break;
				}
				for (k=i+1, l=j+1; k<(ion_cnt1) && l<(ion_cnt2); ++k, ++l) {
					if (spectrum1[k].tag == MAYBE && spectrum2[l].tag == MAYBE && (isless((spectrum1[k].mz - spectrum2[l].mz), TOLERANCE) && isgreater((spectrum1[k].mz - spectrum2[l].mz), -TOLERANCE)))
						make_match(&spectrum1[k], &spectrum2[l], matches1, matches2, &match_cnt);
					else break;
				}
			}
		}
	}
	 
	check_match(spectrum1, ion_cnt1, spectrum2, ion_cnt2, &remaining); /*reset tag as 'none' for those not already matched, then label possible matches 'maybe'*/
	while (remaining) { /*match 'maybes'*/
		for (i=0; i<ion_cnt1; i++) {
                	for (j=0; j<ion_cnt2; j++) { 
				if ((spectrum1[i].tag==MAYBE && spectrum2[j].tag==MAYBE) && (isless((spectrum1[i].mz - spectrum2[j].mz), TOLERANCE) && isgreater((spectrum1[i].mz - spectrum2[j].mz), -TOLERANCE))) {
					make_match(&spectrum1[i], &spectrum2[j], matches1, matches2, &match_cnt);
				}
			}
		}
		check_match(spectrum1, ion_cnt1, spectrum2, ion_cnt2, &remaining);
	}
	*common = match_cnt;
	if (mode != COMMON) { //if want to analyse all ions, add to matches with 0s, if mode = common do nothing
		for (i=0; i<ion_cnt1; i++) { /*add zeros for those with no match*/
			if (spectrum1[i].tag==NONE) {
				memcpy(&matches1[match_cnt], &spectrum1[i], sizeof(ions));
				matches2[match_cnt].mz=0;
				matches2[match_cnt].intensity=0;
				match_cnt++;
			}
		}
		for (i=0; i<ion_cnt2; i++) {
			if (spectrum2[i].tag==NONE) {
				memcpy(&matches2[match_cnt], &spectrum2[i], sizeof(ions));
				matches1[match_cnt].mz=0;
                	        matches1[match_cnt].intensity=0;
				match_cnt++;
			}
		}
	}	
/*print check*/
	fprintf(g, "Matches:\nMatch,MHC m/z,MHC intensity,Synthetic m/z,Synthetic intensity\n");
	for (i=0; i<match_cnt; i++) {
		fprintf(g,"%d,%lf,%lf,%lf,%lf\n",i, matches1[i].mz, matches1[i].intensity, matches2[i].mz, matches2[i].intensity);
	}

	if(((*v1 = (double *)calloc(match_cnt, sizeof(double))) == NULL) || ((*v2 = (double *)calloc(match_cnt, sizeof(double))) == NULL)) {
                printf("Memory allocation error in make_vectors()\n");
                exit(0);
        }
	for (i=0; i<match_cnt; i++) {
		(*v1)[i]=matches1[i].intensity;
		(*v2)[i]=matches2[i].intensity;
	}
	return match_cnt;
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
int sort_by_mz(const void *v1,const void *v2)
{
        ions *m1 = (ions *)v1, *m2 = (ions *)v2;

        if(m1->mz < m2->mz)
                return(-1);

        if(m1->mz > m2->mz)
                return(1);

        return(0);
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void check_match(ions *spectrum1, int ion_cnt1, ions *spectrum2, int ion_cnt2, int *found)
{
	int i, j;
	
	*found=0;
	for (i=0; i<ion_cnt1; i++) {
		if (spectrum1[i].tag != MATCH)
			spectrum1[i].tag = NONE;
	}
	for (i=0; i<ion_cnt2; i++) {
		if (spectrum2[i].tag != MATCH)
			spectrum2[i].tag = NONE;
	}

	for (i=0; i<ion_cnt1; i++) { /*tag all possible matches 'Maybe'*/
		for (j=0; j<ion_cnt2; j++) {
			if (spectrum1[i].tag != MATCH && spectrum2[j].tag != MATCH) {
				if (spectrum1[i].tag == NONE || spectrum2[j].tag == NONE) {
					if ((isless((spectrum1[i].mz - spectrum2[j].mz), TOLERANCE)) && (isgreater((spectrum1[i].mz - spectrum2[j].mz), -TOLERANCE))) {
						spectrum1[i].tag = MAYBE;
						spectrum2[j].tag = MAYBE;
						*found+=1;
					}
				}
			}
		}
	}
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void check_peak(ions *spectrum, int cnt) 
{
	int i, j;

	for (i=0; i<cnt; i++) { /* check spectrum for ions close together and find peak */
		if (spectrum[i].tag == MAYBE) {
			if (i>0 && i<(cnt-1)) {
				if ((isgreater((spectrum[i].intensity/spectrum[(i-1)].intensity), 2.0)) && (isgreater((spectrum[i].intensity/spectrum[(i+1)].intensity), 2.0 ))) {
					spectrum[i].tag = PEAK;
				}
			}
			else if (i==0) {
				if (isgreater((spectrum[i].intensity/spectrum[(i+1)].intensity), 2.0))
					spectrum[i].tag = PEAK;
			}
			else if (i==(cnt-1)) {
				if (isgreater((spectrum[i].intensity/spectrum[(i-1)].intensity), 2.0))
                                        spectrum[i].tag = PEAK;
			} 
		}
	}
	for (i=0; i<cnt; i++) { /*check if a bigger 'maybe' is within tolerance and if so tag that one as the 'peak' instead*/
		if (spectrum[i].tag == PEAK) {
			for (j=(i+1); j<cnt; j++) {
				if (spectrum[j].tag == PEAK && (isless((spectrum[j].mz-spectrum[i].mz), TOLERANCE))) {
					if (isless(spectrum[i].intensity, spectrum[j].intensity)) 
						spectrum[i].tag = MAYBE;
					else spectrum[j].tag = MAYBE;
				}
			}
		}
	}	
	return;
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
void make_match(ions *spectrum1, ions *spectrum2, ions *matches1, ions *matches2, int *match_cnt)
{
	spectrum1[0].tag=MATCH;
	spectrum2[0].tag=MATCH;
	memcpy(&matches1[*match_cnt], spectrum1, sizeof(ions));
	memcpy(&matches2[*match_cnt], spectrum2, sizeof(ions));
	*match_cnt+=1;
}
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
double dotproduct(double *v1, double *v2, int cnt)
{
	int i;
	double sum = 0.0;

	for (i=0; i<cnt; ++i) {
		sum += v1[i] * v2[i];
	}
	return sum;
}
/*****************************************************************************************************************************/




