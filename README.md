# immunopeptidogenomics
## Tools for harnessing RNA-seq data to discover cryptic peptides in the immunopeptidome by mass spectrometry
This repository contains source code for compiling the novel software needed in our immunopeptidogenomics workflow, which produces sample-specific, 3-frame translated transcriptomics-based protein databases for identifying peptides from tandem mass spectrometry data, and further tools to help decipher the search results.

Most code is in C and was produced and compiled on Linux (Ubuntu 16.04). It should also be possible to compile for Windows but this has not been tested. Further details regarding each tool and instructions for compilation and usage on Linux follow.

More details and usage examples can be found in our paper (manuscript under review as of 18/12/2020)

### alt_liftover
This tool performs 'liftover' for gtf files which refer to a standard reference genome (e.g.  GRCh38), so that the output gtf file contains coordinates referencing an alternate genome (produced by GATK’s FastaAlternateReferenceMaker from the same reference genome and a variant vcf file).

##### NOTE: USE A CURATED VCF FILE (e.g. produced by curate_vcf.c) so that you can be sure what FastaAlternateReferenceMaker has done!

Compilation example:

```
cc alt_liftover.c –o alt_liftover.o
```
Usage example:
```
./alt_liftover.o –v vcf_file.vcf –g gtf_file.gtf –f _output &> my.log
```
|Required Input:||
|---|---|
|-v|curated vcf file (previously used to make alternate genome)|
|-g|gtf file (sorted so that entries for each chromosome are grouped together)|

|Optional:||
---|---
-f|Use to specify a suffix for the output gtf file (default is ‘_alt’)
-h|print help

**Output**: liftover gtf file 

### curate_vcf
This tool curates a vcf file in preparation to provide an unambiguous input.vcf containing simple substitutions and indels for GATK’s FastaAlternateReferenceMaker, so the user knows exactly which variants FastaAlternateReferenceMaker will include in the alternate genome. The intention is to run the program twice, with and without –d option.

Compilation example:
```
cc curate_vcf.c –o curate_vcf.o
```
Usage example:
```
./curate_vcf.o [-d] input.vcf 
```

**Required input:** vcf file. 

|Optional:||
---|---
-d|use this to 'deprioritise deletions' i.e. remove from vcf any deletion mutations which would mask downstream mutation sites

**Output**: curated vcf file named (input_name)_unmasked.vcf or (input_name)_indel.vcf, depending on whether or not the –d option is selected.

### msDot
This tool judges the similarity in fragmentation pattern between pairs of spectra by calculating normalised dot products. It compares pairs of dta files following the naming convention convention x_mhc_peptide.dta and x_synth_peptide.dta, where x is an identifier such as the presumed peptide sequence. The program reads in the ions, cleans spectra to remove peaks with intensity below 0.01 * max intensity, and makes vectors with equal numbers of values in each array (i.e. same number of peaks), matching peaks within 0.1 Da between spectra. Then for vectors a and b,

Normalised dot product = a∙b / (|a||b|) where |a| = √(a∙a) and |b| = √(b∙b)

This yields a value between 0 and 1 where closer to 1 means more similar.

Compilation example:
```
cc msDot.c -o msDot.o –lm
```
Usage example:
```
./msDot.o input_folder 
```
**Input:** 	Path to directory containing pairs of dta files named to follow convention IDENTIFIER_mhc_peptide.dta and IDENTIFIER_synth_peptide.dta
|Output files (in same directory):||
---|---
dotproducts_full.csv|extensive information regarding dot product calculation including peak matching between spectra 
dotproducts_only.csv|simple table of results for each pair, including dot product and the ratio of the (number of ions common to both spectra)/(least number of ions per spectrum in pair)

### origins
This tool interrogates your transcriptomes, standard protein database and Ensembl to provide information about the possible origins of unconventional sequences. 

Compilation example:

Origins has dependencies on curl and cJSON. You will need to install/download these. For example, on Ubuntu curl can be installed using `sudo apt install curl`

cJSON can be sourced from https://github.com/DaveGamble/cJSON; clone the repository/download the cJSON.h and cJSON.c source files. Next, open origins.c and alter line 9 to include the filepath to cJSON.h on your system. Then compile origins including cJSON.c, e.g.
```	 
cc origins.c -o origins.o `curl-config --libs` /path/to/cJSON.c
```
Usage example:
```
./origins.o -p peptides.txt \
-r cuffcompare_output.tracking \
-d uniprot.fasta \
-n no_variant_transcriptome.fa \
[-a indel_transcriptome.fa,unmasked_transcriptome.fa\
 	-v indel.vcf,unmasked.vcf]
```
|Required Input:||
---|---
-p|txt file listing peptide sequences to check
-r|Cuffcompare output tracking file 
-d|standard protein database (e.g. Uniprot)
-n|‘normal’ transcriptome i.e. without variants added (nucleotide sequences in fasta format) 

|Optional (see note):||
---|---
-a|comma-separated list of ‘alternate’ transcriptomes i.e. without variants added (nucleotide sequences in fasta format)
-v|comma-separated list of the vcf files used to produce the corresponding transcriptomes
-h|print help

Note: -a and –v are optional, but if used the files should correspond; e.g. `–a transcriptome1.fa,transcriptome2.fa –v v1.vcf,v2.vcf`

|Output:||
---|---
x_origins_rna.csv|file of peptide sequences with possible transcripts, and for each transcript, coordinates and metadata
x_origins_prot.csv|file of only the input peptide sequences which might have a conventional source from the (-d) protein database, with protein details 
x_origins_discard.txt|file of sequences found neither in the translated transcriptomes nor the protein database, hence artificial sequences to be discarded 

where x is the title of the peptide list input file

### revert_headers
FastaAlternateReferenceMaker changes the chromosome names; this tool changes them back to match the reference genome. It outputs the revised genome as tmp.fasta so that the user can manually check in the log output that it has matched chromosomes correctly before overwriting the original version. (This should be fine if you use the same reference genome used with FastaAlternateReferenceMaker, but GATK might change their program so it’s worth checking.)

Compilation example:
```
cc revert_header.c –o revert_headers.o &> my.log
```
Usage example:
```
./revert_headers.o reference_genome.fa alternate_genome_to_fix.fa
```
**Required input IN ORDER:** reference_genome.fa THEN alternate_genome_to_fix.fa

**Output:** tmp.fasta

Also the stdout will list the changes made e.g. “Was x now y”

### squish
This tool is designed to merge the protein databases produced by running triple_translate on various transcriptomes, and minimise redundancy by removing duplicates and sequences wholly contained within another sequence. 

Compilation example:
```
cc squish.c -lm -lpthread -o squish.o
```
Usage example:
```
./squish.o -d prot_db1.fasta -d prot_db2.fasta [-t 8] [-o another_name.fasta] &> squish.log
```
|Required input:||
---|---
-d|specify one protein fasta database (this option should be repeated to add more databases; e.g. `–d database1.fasta –d database2.fasta`)

|Optional:||
---|---
-t|Specify number of threads to use (default = 1)
-o|Specify name of output database file (default = output.fasta)
-h|print usage help

### triple_translate
This tool takes a transcriptome file of nucleotide sequences in fasta format and performs 3 frame translation. Optionally, the cuffcompare tracking file produced during transcriptome assembly can be input to include associated gene information in the protein sequence headers. triple_translate prints any sequences > 7 aa in length to a new protein database in fasta format.

Compilation example:
```
cc triple_translate.c –o triple_translate.o
```
Usage example:
```
./triple_translate.o [-c cuffcompare_output.tracking] transcriptome.fa
```
**Required input:** specify one transcriptome sequences .fa/.fasta file for translation

|Optional:||
---|---
-c|Cuffcompare output tracking file
-h|print usage help

**Output:**	protein fasta file with name x_3translate.fasta where x is the input file name

### db_compare.R
This tool helps compare ‘DB search psm.csv’ results from a PEAKS search on the same data against two databases (labelled standard and cryptic), outputting various graphs and lists.

**Dependencies:** R packages VennDiagram, ggplot2, dplyr, hrbrthemes, reshape2, optparse

Usage example:
```
Rscript db_compare.R –c cryptic_psms.csv –d 15.0 –n normal_psms.csv –m 14.0 [-p output_prefix] [-j origins_discard.txt]
```
|Options:||
---|---
-c CHARACTER, --cryptic=CHARACTER|cryptic PEAKS results filename
-n CHARACTER, --normal=CHARACTER|normal (e.g. uniprot) PEAKS results filename
-j CHARACTER, --junction=CHARACTER|Optional: artificial junction peptide list for discard (txt)
-p CHARACTER, --prefix=CHARACTER|prefix for all output files (e.g. myCellsRep1)
-d NUMBER, --cryptic_threshold=NUMBER|threshold score for cryptic search results
-m NUMBER, --norm_threshold=NUMBER|threshold score for normal (uniprot) search results
-h, --help|Show this help message and exit
