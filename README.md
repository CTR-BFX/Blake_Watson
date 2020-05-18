# Genetic disruption of folate metabolism causes epigenetic instability in sperm and distinguishes HIRA as a biomarker of maternal epigenetic inheritance
**Georgina E.T. Blake<sup>1,2,†</sup>, Xiaohui Zhao<sup>1,2</sup>, Hong wa Yung<sup>1,2</sup>,Graham J. Burton<sup>1,2</sup>,Anne C.Ferguson-Smith<sup>2,3</sup>,
Russell Hamilton<sup>2,3</sup>, and Erica D. Watson<sup>1,2*</sup>** <br>
<sup>1</sup> Department of Physiology, Development and Neuroscience, University of Cambridge, Cambridge, UK <br>
<sup>2</sup> Centre for Trophoblast Research, University of Cambridge, Cambridge UK<br>
<sup>3</sup> Department of Genetics, University of Cambridge, Cambridge, UK<br>
<sup>†</sup> Current address: College of Medicine and Health, University of Exeter Medical School, Exeter, UK<br>
<sup>*</sup> corresponding author: edw23@cam.ac.uk<br>

### Publication ###

Georgina E.T. Blake, Xiaohui Zhao, Hong wa Yung, Graham J. Burton, Anne C.Ferguson-Smith, Russell S. Hamilton and Erica D. Watson (2020) Genetic disruption of folate metabolism causes epigenetic instability in sperm and distinguishes HIRA as a biomarker of maternal epigenetic inheritance.

Code Release to accompany paper:


## Whole Genome Sequencing Analysis (WGS) <br>
 - [x] Samples: C57BI/6J (Control-1, Control-2), Mtrr<sup>{gt/gt}</sup> (800-1, 800-2, 800-3, 800-4, 800-5 and 800-6)
 - [x] SNPs (single nucleotide polymorphisms) and SVs (structure variants) calling for C57BI/6J control embryos and Mtrr<sup>{gt/gt}</sup> embryos

### Sequencing Data Quality Control (FastQC, v0.11.5), Adapter trimming (Trim_galore, v0.6.4) and report summary (MultiQC, v.1.4);


### Sequencing Aligned to the C57Bl/6J reference genome(GRCm38, mm10) using Bowtie2(v2.3.4) with default parameters;

Alignments processed using ClusterFlow [[GitHub](https://github.com/ewels/clusterflow)] [[DOI](http://dx.doi.org/10.12688/f1000research.10335.2)]

Example command line as used in clusterflow:

    READ1="160326_X169_FCHMVV5CCXX_L4_8521603000353_1_val_1.fq.gz"
    READ2="160326_X169_FCHMVV5CCXX_L4_8521603000353_2_val_2.fq.gz"

    bowtie2 -p 4 -t --phred33-quals \
            -x Genomes/Mus_musculus/GRCm38/GRCm38 \
            -1 ${READ1} -2 ${READ2} | \
            samtools view -bS - > ${READ1/.fq.gz/}_bowtie2.bam

> Duplicates were marked using Picard tools (v2.24.0);

    READ1="160326_X169_FCHMVV5CCXX_L4_8521603000353_1_val_1.fq.gz"
    READ2="160326_X169_FCHMVV5CCXX_L4_8521603000353_2_val_2.fq.gz"

    java -Xmx2g -jar picard-tools-2.2.4/picard.jar MarkDuplicates \
         INPUT=${READ1/.fq.gz/}_bowtie2.bam \
         OUTPUT=${READ1/.fq.gz/}_bowtie2.bam_MarkDups.bam \
         ASSUME_SORTED=false \
         REMOVE_DUPLICATES=false \
         METRICS_FILE=1${READ1/.fq.gz/}_bowtie2.bam_picardDupMetrics.txt \
        VALIDATION_STRINGENCY=LENIENT

### Mtrr region masking;

    PROJECT="CTR_edw23_0001"
    BEDFILE="chr13-mtrr-mask-region_20Mb.bed"
    BAMSTEM="srtd_mrgd_srtd_dedup_chr13-mtrr-20Mb-masked"

    declare -a samples=("Control-1" "Control-2" "800-1" "800-2" "800-3" "800-4" "800-5" "800-6")

    for i in "${samples[@]}"
        do
          echo "sample = $i  "
          cd ${i}
          samtools view ${PROJECT}_${i}.srtd_mrgd_srtd_dedup.bam \
                   -b -h \
                   -o ${PROJECT}_${i}.${BAMSTEM}.bam \
                   -U ${PROJECT}_${i}.${BAMSTEM}.bam \
                   -L ${BEDFILE}
        done

### Run Manta Structural Variant Caller

Manta (v0.29.6) [[GitHub](https://github.com/Illumina/manta)] [[DOI](https://doi.org/10.1093/bioinformatics/btv710)]

> Set up some variables for running Manta

    PROJECT="CTR_edw23_0001"
    GENOME="Genomes/Mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.chromosome.all.fasta"
    BAMSTEM="srtd_mrgd_srtd_dedup_chr13-mtrr-20Mb-masked"
    MANTAWORKDIR="MantaWorkflow_20Mb_Masked"

    declare -a samples=("Control-1" "Control-2" "800-1" "800-2" "800-3" "800-4" "800-5" "800-6")


##### Run Manta Prepare

    for i in "${samples[@]}"
      do
        echo "sample = $i  "
        mkdir -p ${i}/${MANTAWORKDIR}
        configManta.py --bam ${PROJECT}_${SAMPLE}.${BAMSTEM}.bam --referenceFasta ${GENOME} --runDir ${i}/${MANTAWORKDIR}
     done

##### Run Manta workflow

    for i in "${samples[@]}"
      do
        echo "sample = $i  "
        cd ${i}/${MANTAWORKDIR}
        runWorkflow.py --mode=local --jobs=10 --memGb=25
      done

### Remapped to the mm10 reference mouse genome using BWA (v480 0.7.15-r1144- dirty);

Alignments processed using ClusterFlow [[GitHub](https://github.com/ewels/clusterflow)] [[DOI](http://dx.doi.org/10.12688/f1000research.10335.2)]

Example command line as used in clusterflow:

    READ1="160326_X169_FCHMVV5CCXX_L4_8521603000353_1_val_1.fq.gz"
    READ2="160326_X169_FCHMVV5CCXX_L4_8521603000353_2_val_2.fq.gz"

    bwa mem -t 12 \
            Genomes/Mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.chromosome.all.fa.gz \
            ${READ1} ${READ2} |\
            samtools view -bS - > ${READ1/.fq.gz/}_GRCm38_bwa.bam

> Duplicates were marked using Picard tools (v2.24.0)

    java -Xmx2g -jar picard-tools-2.2.4/picard.jar MarkDuplicates \
         INPUT=${READ1/.fq.gz/}_GRCm38_bwa.bam \
         OUTPUT=${READ1/.fq.gz/}_GRCm38_bwa.bam_MarkDups.bam \
         ASSUME_SORTED=false REMOVE_DUPLICATES=false \
         METRICS_FILE=1${READ1/.fq.gz/}_GRCm38_bwa.bam_GRCm38_bwa_srtd.bam_picardDupMetrics.txt \
         VALIDATION_STRINGENCY=LENIENT        

### Reads were locally realigned and SNPs and short indels identified using GenomeAnalysisTK (GATK, v3.7) <br>

> Set up some variables

    REFERENCE_FASTA="Genomes/Mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.chromosome.all.fasta"
    KNOWNSITES_1="KNOWNSNPs/dbSNP_vcf_chr_ALL.vcf.gz"
    KNOWNSITES_2="KNOWNSNPs/129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz"
    GATKLOCN="GenomeAnalysisTK-3.8-0"
    PICARDLOCN="picard-tools-2.2.4/"
    TMPDIR=$(pwd)

##### Recalibrate Base Quality Scores

> Example BAM file

    INPUTBAM="Control-1.GRCm38_bwa_srtd.bam_MarkDups.RG.bam"
    INPUTBAMPRERG=${INPUTBAM/MarkDups.RG.bam/MarkDups.bam}

    SAMPLE=${INPUTBAM/.GRCm38_bwa_srtd.bam_MarkDups.bam/}
    ID="ID:"${SAMPLE}
    SM="SM:"${SAMPLE}
    LB="LB:"${SAMPLE}

> AddOrReplaceReadGroups

    java -Djava.io.tmpdir=${TMPDIR} -Xmx5g -jar ${PICARDLOCN}/picard.jar AddOrReplaceReadGroups \
         I=${INPUTBAMPRERG} \
         O=${INPUTBAM} \
         RGID=${SAMPLE} \
         RGLB=${SAMPLE} \
         RGPL=illumina \
         RGPU=${SAMPLE} \
         RGSM=${SAMPLE}

>  Index BAM file

    java -jar ${PICARDLOCN}/picard.jar BuildBamIndex \
         I=${INPUTBAM} \
         TMPDIR=${TMPDIR}

> Analyze patterns of covariation in the sequence dataset

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T BaseRecalibrator \
         -nct 8 \
         -R ${REFERENCE_FASTA} \
         -I ${INPUTBAM} \
         -knownSites ${KNOWNSITES_1} \
         -knownSites ${KNOWNSITES_2} \
         -o ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recal_data.table}


> Do a second pass to analyze covariation remaining after recalibration

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T BaseRecalibrator \
         -nct 8 \
         -R ${REFERENCE_FASTA} \
         -I ${INPUTBAM} \
         -knownSites ${KNOWNSITES_1} \
         -knownSites ${KNOWNSITES_2} \
         -BQSR ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recal_data.table} \
         -o ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.post_recal_data.table}

> Generate before/after plots

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T AnalyzeCovariates \
         -l DEBUG \
         -R ${REFERENCE_FASTA} \
         -before ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recal_data.table} \
         -after ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.post_recal_data.table} \
         -plots ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recalibration_plots.pdf}

> Apply the recalibration to your sequence data

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T PrintReads \
         -nct 8 \
         -R ${REFERENCE_FASTA} \
         -I ${INPUTBAM} \
         -BQSR ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recal_data.table} \
         -o ${INPUTBAM/MarkDups.RG.bam/MarkDups.RG.recal.bam}

##### Call Variants

    INPUTBAM="Control-1.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.bam"

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R ${REFERENCE_FASTA} \
         -I ${INPUTBAM} \
         --variant_index_type LINEAR \
         --variant_index_parameter 128000 \
         -ERC GVCF \
         -o ${INPUTBAM/recal.bam/recal.raw_variants.vcf}

##### Joint Genotyping across all samples

    java -Djava.io.tmpdir=${TMPDIR} -jar ${GATKLOCN}/GenomeAnalysisTK.jar \
         -T GenotypeGVCFs \
         -R ${REFERENCE_FASTA} \
         --variant 800-1/800-1.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant 800-2/800-2.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant 800-3/800-3.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant 800-4/800-4.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant 800-5/800-5.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant 800-6/800-6.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant Control-1/Control-1.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         --variant Control-2/Control-2.GRCm38_bwa_srtd.bam_MarkDups.RG.recal.raw_variants.vcf \
         -o joint_genotyping_output.vcf


### SNPs filter using vcftools (vcf-annotate), [Merge_Filter_Snp_Annotation.sh](https://github.com/xz289/Blake_Watson/blob/master/Scripts/Merge_Filter_Snp_Annotation.sh) <br>
       - [x] 1-StrandBias [0.0001] ---Min P-value for strand bias
       - [x] 2-BaseQualBias [0.0002] ---Min P-value for BaseQ bias
       - [x] 3-MapQualBias [0.00001] ---Min P-value for mapQ bias
       - [x] 4-EndDistBias [0.0001] ---Min P-value for end distance bias
       - [x] a-MinAB [6] ---Min number of alternate bases
       - [x] d-MinDP [14] ---Min read depth
       - [x] D-MaxDP [100] ---Max read depth
       - [x] q-MinMQ [25] ---Min RMS(Root Mean Square) mapping quality for SNPs
       - [x] Q-Qual [40] ---Min value of the QUAL field
       - [x] W-GapWin [3] ---Window size for filtering adjacent gaps
       - [x] w-SnpGap [5] ---SNP within 5bp around a gap to be filtered
       - [x] H-HWE [0.0001] ---Min P-value for HWE(Hardy–Weinberg equilibrium) and (F\<0)
       - [x] v-VDB [0] ---Min Variant Distance bias

### Extra filtering the repeatitive and defind Heterozygous SNPs. R script [Exclued_SNPs_Homo_Dinu_SegDup_Het.R](Scripts/Exclued_SNPs_Homo_Dinu_SegDup_Het.R).
      Reference: Oey, H., Isbel, L., Hickey, P., Ebaid, B. & Whitelaw, E. Genetic and epigenetic variation among inbred mouse littermates: identification of inter-individual differentially methylated regions. Epigenetics Chromatin 8, 54 (2015). <br>
       - [x] Simple repeats periodicity < 9bp;
       - [x] Homopolymer repeats >8 bp;
       - [x] dinucleotide repeats >14 bp;
       - [x] low mapping quality (<40);
       - [x] Heterozygous variant call counts > 3 within 10kb window, overlap with segment duplication or repeats;


### Corresponding Figures (Figure 1(A), Supplementary Figure 3 (C), (D)) generation Code is <br>
R script [SV_SNPs_Fig1A_SFig3C_3D.R](Scripts/SV_SNPs_Fig1A_SFig3C_3D.R)

## MeDIP Data Analysis <br>

##### Make genome interval files

[make_genome_intervals.pl](MeDIPs/make_genome_intervals.pl)
creates intervals of set sizes for a genome of interest. Requires a fai genome index file for the genome created with e.g. `samtools faidx genome.fasta`. Output is a BED format file

    perl make_genome_intervals.pl 5000 \
    Mus_musculus.GRCm38.dna.chromosome.all.fa.fai  | \
    sort -k 1,1n -k2,2n > output.bed

##### Make Bedtools Coverage Histogram

[makeBedtoolsCoverageHist.sh](MeDIPs/makeBedtoolsCoverageHist.sh) takes the sample bam (ending in .srtd.bam) files in a directory and runs bedtools coverage for the defined genomic interval files (hard coded, edit file to change)

    ./makeBedtoolsCoverageHist.sh

The resulting bed file is the output from `bedtools coverage` with the following columns"

 - The number of features in the bam file that overlapped (by at least one base pair) the genome intervals.
 - The number of bases in genome intervals that had non-zero coverage from features in bam file.
 - The length of the entry in genome intervals.
 - The fraction of bases in genome intervals that had non-zero coverage from features in the bam file.

##### Combine Interval Tables

 [combineTables.pl](MeDIPs/combineTables.pl) and [combineTables_command.sh](MeDIPs/combineTables_command.sh) take all the individual sample coverage bed files e.g. `C57_36.chr.srtd.25Kb.bed` in a directory and makes a single combined table. The interval size in Kb is given on the command line.

     perl combineTables.pl 25 > combineTables_command.sh

 The above command produces a bash script of the normalised coverage per genomic window doe each of the samples

     ./combineTables_command.sh

 Running the shell script produces the final table for analysis in R

 Output is a list of commands to be run on the command line. This step is only required if the samples to be analysed changes. i.e. samples removed or added.

 A pre-made version of the commands available in `combineTables_command.sh`

     ./combineTables_command.sh

 Output is a table `CTR_edw23_0003_MeDIPS_Region_Coverage_Combo_${WindowSize}Kb.table`

 ##### Create PCA Plots

 [CTR_edw23_0003_MeDIPS.R](MeDIPs/CTR_edw23_0003_MeDIPS.R) script for calculating PCA plots for genomic interval coverage

 ##### Template Script
 [GB_Template.R](MeDIPs/GB_Template.R)
 Template R script for performing MEDIPs analysis

 #####  Convert BAM Chr labels
 [convertBAM_chr.sh](MeDIPs/convertBAM_chr.sh) a simple script to convert chromosome nomenclature in BAM files from 1 to chr1. Utilises samtools view and reheader


## ChIP-Seq and ATAC-Seq Analysis <br>

### Step1: Data availability(Supplementary files 2), relating published reference:

#### Mouse: Epiblast (Epi) and extraembryonic ectoderm (ExE) from E6.5,  ChIPSeq (GSE84236, mm9, THSS)
* Zachary D. smith, Jiantao shi, Hongcang Gu, Julie Donaghey, Kendell Clement, Davide Cacchiarelli, Andreas Gnirke, Franziska michor & Alexander meissner [[DOI](doi: 10.1038/nature23891)]

#### Mouse: Sperm ChIPSeq  (GSE79230, mm9)
* Jung YH, Sauria MEG, Lyu X, Cheema MS, Ausio J, Taylor J, Corces VG.[DOI](doi: 10.1016/j.celrep.2017.01.034.)]

#### Mouse: trophoblast stem cells (TSCs) and embryonic stem cells (ESCs) ChIPSeq (Supplementary Files 2, mm10)
* Stefan Schoenfelder, Borbala Mifsud, Claire E. Senner, Christopher D. Todd,Stephanie Chrysanthou, Elodie Darbo, Myriam Hemberger & Miguel R. Branco[[DOI](DOI: 10.1038/s41467-018-06666-4)]


### Step2: Remove 20Mb region surrounding the Mtrr gene (Chr13:58060780-80060780)
DMRs total number reduced from 893 to 459, (Mtrr<sup>{+/gt}</sup>, Mtrr<sup>{+/+}</sup> and Mtrr<sup>{gt/gt}</sup>)
> Table: Number of summary DMRs

| Genotype   |  fileName | no. of DMRs|  
| ----- | --- |--- |
| Mtrr<sup>{+/gt}</sup>  |  [merge_C57_Vs_+gt_0.01_tab1_dmr_xy_29418_allcol.bed](MeDIPs/merge_C57_Vs_+gt_0.01_tab1_dmr_xy_29418_allcol.bed) | 203 |
| Mtrr<sup>{+/+}</sup>   |  [merge_C57_Vs_++_0.01_xy_6518_tab1_dmr_allcol.bed](MeDIPs/merge_C57_Vs_++_0.01_xy_6518_tab1_dmr_allcol.bed) | 91 |
| Mtrr<sup>{gt/gt}</sup> |  [merge_C57_Vs_gtgt_0.01_xy_allcol.bed](MeDIPs/merge_C57_Vs_gtgt_0.01_xy_allcol.bed) | 599|
| Mtrr.remove |  [AllDMRs_maskMtrr20Mb_N459.bed](MeDIPs/AllDMRs_maskMtrr20Mb_N459.bed) | 459|

### Step 3: hgLiftOver the mm9 ChIP/ATAC bw file to mm10, and clip the extra bps for MACS output, then compute profiles.
| Script   |  fileName |  
| ----- | --- |
| bash  |  [BedClip_MACS_overlap.sh](Scripts/BedClip_MACS_overlap.sh) |
| bash  |  [HgLiftOver_RandomDMRs_ProfilePlot.sh](Scripts/HgLiftOver_RandomDMRs_ProfilePlot.sh) |

> Profiles Plots (see Supplementary Fig 7)

### Step 4: Selected candidates genome viewer using IGV(v2.5.2). <br>

| Gene  | PNG |
| ----- | --- |
| Hira | [Hira_refined_IGV_snapshot.png](Figures/Hira_refined_IGV_snapshot.png) |
| Cwc27 | [Cwc27_refined_IGV_snapshot.png](Figures/Cwc27_refined_IGV_snapshot.png) |
| Tshz3  | [Tshz3_refined_IGV_snapshot.png](Figures/Tshz3_refined_IGV_snapshot.png) |

## Data AarrayExpress Link
| DataType  | Accession | Link |
| ----- | --- | ---|
| Mtrr MeDIP Data | E-MTAB-8533 | [E-MTAB-8533](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8533/) |
| Mtrr WGS Data | [E-MTAB-8513]|E-MTAB-8513](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8513) |



````

## Contact

Contact Xiaohui Zhao (xz289 -at- cam.ac.uk) and Russell S. Hamilton (rsh46 -at- cam.ac.uk)
