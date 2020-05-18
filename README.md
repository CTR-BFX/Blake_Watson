# Genetic disruption of folate metabolism causes epigenetic instability in sperm and distinguishes HIRA as a biomarker of maternal epigenetic inheritance
**Georgina E.T. Blake<sup>1,2,†</sup>, Xiaohui Zhao<sup>1,2</sup>, Hong wa Yung<sup>1,2</sup>,Graham J. Burton<sup>1,2</sup>,Anne C.Ferguson-Smith<sup>2,3</sup>,
Russell Hamilton<sup>2,3</sup>, and Erica D. Watson<sup>1,2*</sup>**
<sup>1</sup> Department of Physiology, Development and Neuroscience, University of Cambridge, Cambridge, UK <br>
<sup>2</sup> Centre for Trophoblast Research, University of Cambridge, Cambridge UK<br>
<sup>3</sup> Department of Genetics, University of Cambridge, Cambridge, UK<br>
<sup>†</sup> Current address: College of Medicine and Health, University of Exeter Medical School, Exeter, UK<br>
<sup>*</sup> corresponding author: edw23@cam.ac.uk<br>

### Publication ###

Georgina E.T. Blake, Xiaohui Zhao, Hong wa Yung, Graham J. Burton, Anne C.Ferguson-Smith, Russell S. Hamilton and Erica D. Watson (2020) Genetic disruption of folate metabolism causes epigenetic instability in sperm and distinguishes HIRA as a biomarker of maternal epigenetic inheritance. [[xx]](http://www.xx) [[DOI]](https://doi.org/xx)

Code Release to accompany paper: [![DOI](xx)](xx)


## Whole Genome Sequencing Analysis (WGS) <br>
 C57BI/6J (Control-1, Control-2), Mtrr^{gt/gt} (800-1, 800-2, 800-3, 800-4, 800-5 and 800-6)
### Part I: Whole Genome Sequencing (WGS) SNPs (singlenucleotide polymorphisms) calling for C57BI/6J control embryos and Mtrr^{gt/gt} embryos
#### Step 1: Sequencing Data Quality Control (FastQC, v0.11.5), Adapter trimming (Trim_galore, v0.6.4) and report summary (MultiQC, v.1.4); ;
#### Step 2: Sequencing Aligned to the C57Bl/6J reference genome(GRCm38, mm10) using Bowtie2(v2.3.4) with default parameters;
#### Step 3: Duplicates were marked using Picard tools (v2.9.0);
#### Step 4: Remapped to the mm10 reference mouse genome using BWA (v480 0.7.15-r1144- dirty);
#### Step 5: Reads were locally realigned and SNPs and short indels identified using GenomeAnalysisTK (GATK, v3.7) <br>
     - [x] MarkDuplicates
     - [x] HaplotypeCaller
     - [x] SelectVariants
     - [x] CombineVariants
     - [x] VariantRecalibrator
     - [x] ApplyRecalibration
#### Step 6: SNPs filter using vcftools (vcf-annotate,[a link][a link] https://github.com/xz289/Blake_Watson/blob/master/Scripts/Merge_Filter_Snp_Annotation.sh ) <br>
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
#### Filtering the repeatitive and defind Heterozygous SNPs. [a link] https://github.com/xz289/Blake_Watson/blob/master/Scripts/Exclued_SNPs_Homo_Dinu_SegDup_Het.R <br>
      Reference: Oey, H., Isbel, L., Hickey, P., Ebaid, B. & Whitelaw, E. Genetic and epigenetic variation among inbred mouse littermates: identification of inter-individual differentially methylated regions. Epigenetics Chromatin 8, 54 (2015). <br>
       - [x] Simple repeats periodicity < 9bp;
       - [x] Homopolymer repeats >8 bp;
       - [x] dinucleotide repeats >14 bp;
       - [x] low mapping quality (<40);
       - [x] Heterozygous variant call counts > 3 within 10kb window, overlap with segment duplication or repeats;
### Part II: Whole Genome Sequencing (WGS) SV (structure variant) calling for C57BI/6J control embryos and Mtrr^{gt/gt} embryos


## MeDIP Data Analysis <br>

## ChIP-Seq and ATAC-Seq Analysis



````

## Contact

Contact Xiaohui Zhao (xz289 -at- cam.ac.uk) and Russell S. Hamilton (rsh46 -at- cam.ac.uk)
