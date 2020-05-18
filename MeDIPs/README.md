# CTR_edw23_0003 #


**Russell S. Hamilton** and **Gina Blake**

*Centre for Trophoblast Research, Department of Physiology, Development and Neuroscience, University of Cambridge, Downing Site, Cambridge, CB2 3DY, UK*

GPLv3 ::: September 2017


## Scripts ##
These are not intended as generic analysis tools so require editing / configuring as indicated.

### `make_genome_intervals.pl` ###
Creates intervals of set sizes for a genome of interest. Requires a fai genome index file for the genome created with e.g. `samtools faidx genome.fasta`

    perl make_genome_intervals.pl 5000 Mus_musculus.GRCm38.dna.chromosome.all.fa.fai  | sort -k 1,1n -k2,2n > output.bed

Output is a BED format file   

### `makeBedtoolsCoverageHist.sh` ###
Takes the sample bam (ending in .srtd.bam) files in a directory and runs bedtools coverage for the defined genomic interval files (hard coded, edit file to change)

    ./makeBedtoolsCoverageHist.sh

The resulting bed file is the output from `bedtools coverage` with the following columns"

 - The number of features in the bam file that overlapped (by at least one base pair) the genome intervals.
 - The number of bases in genome intervals that had non-zero coverage from features in bam file.
 - The length of the entry in genome intervals.
 - The fraction of bases in genome intervals that had non-zero coverage from features in the bam file.


### `combineTables.pl` and `combineTables_command.sh` ###
Takes all the individual sample coverage bed files e.g. `C57_36.chr.srtd.25Kb.bed` in a directory and makes a single combined table. The interval size in Kb is given on the command line.

    perl combineTables.pl 25 > combineTables_command.sh

The above command produces a bash script of the normalised coverage per genomic window doe each of the samples

    ./combineTables_command.sh

Running the shell script produces the final table for analysis in R

Output is a list of commands to be run on the command line. This step is only required if the samples to be analysed changes. i.e. samples removed or added.

A pre-made version of the commands available in `combineTables_command.sh`

    ./combineTables_command.sh

Output is a table `CTR_edw23_0003_MeDIPS_Region_Coverage_Combo_${WindowSize}Kb.table`

### `CTR_edw23_0003_MeDIPS.R` ###
Script for calculating PCA plots for genomic interval coverage

### `GB_Template.R` ###
Template R script for performing MEDIPs analysis

### `convertBAM_chr.sh` ###
Simple script to convert chromosome nomenclature in BAM files from 1 to chr1. Utilises samtools view and reheader
