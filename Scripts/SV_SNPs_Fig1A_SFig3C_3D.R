#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# WGS SNPs, SVs annotation, and Figures production in the paper  
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/Blake_Watson
#
# 
# Analysis Performed by Xiaohui Zhao
# CTR Bioinformatics Facility
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+---------------Basic settings and libraries calling----------------------------------------+")

suppressPackageStartupMessages({
library(VariantAnnotation)  ## (v1.28.13) mm10 Ref
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(StructuralVariantAnnotation)
library(segmentSeq)
library(ggplot2)
library(reshape)
library(vcfR)
library(UpSetR)
library(dplyr)
library(cowplot)
library(GenomicRanges)
})

setwd("/storage/CTR-Projects/CTR_edw23/CTR_edw23_0001")
Project <- "CTR_edw23_0001"
vcf.names <- c("800-1", "800-2", "800-3", "800-4", "800-5", "800-6", "Control-1", "Control-2")
vcf.files <- paste0(Project, "_", vcf.names,"_Manta_20MbMasked.diploidSV.vcf")

message("+--------------------VariantAnnotation calling----------------------------------------------+")

Freq.table <- list()
for(a in 1:8){
  vcf <- readVcf(paste0("./Manta_20Mb_Masked_Results/", vcf.files[a]))
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  rd <- rowRanges(vcf)
  seqlevels(rd)[20] <- "M"
  seqlevels(rd) <- paste0("chr",seqlevels(rd)) 
  loc <- locateVariants(rd, txdb, AllVariants())
  Freq.table[[a]] <- as.data.frame(table(loc$LOCATION));
  Freq.table
}

Freq.tableNew <- cbind(Freq.table[[1]], Freq.table[[2]][,2],Freq.table[[3]][,2], Freq.table[[4]][,2],
                       Freq.table[[5]][,2],Freq.table[[6]][,2], Freq.table[[7]][,2],
                       Freq.table[[8]][,2])
colnames(Freq.tableNew) <- c("Annotation", vcf.names)
Freq.total <- colSums(Freq.tableNew[,-1])
rownames(Freq.tableNew) <- Freq.tableNew[,1]
Freq.tableNew <- Freq.tableNew[,-1]
Freq.tableNew <- t(Freq.tableNew)
Freq.tableProp <- Freq.tableNew/Freq.total*100

write.csv(Freq.tableProp, file=paste0("./Manta_20Mb_Masked_Results/", Project, "-SVs_R.VariantAnnotation_Summary_Table_Proportion_07_05_19.csv"))
write.csv(Freq.tableNew, file=paste0("./Manta_20Mb_Masked_Results/", Project, "-SVs_R.VariantAnnotation_Summary_Table_Number_07_05_19.csv"))



message("+--- R calling vcftools to split merged snp annotation file into subsets by samples --------+")

system("vcf-subset -c  recalibrated_vcf_filter1_segDup_simpRL9_hete3_RefH8D14_GRCm3884_ud0.ann.vcf.gz | bgzip -c > out.vcf.gz")

message("+------R calling ANNOVAR to run annotation --------------------------------------+")

GRDataFrame <- list()
GRDataManta.tab <- matrix(0, nrow = 5, ncol = 1)
rownames(GRDataManta.tab) <- c("BND", "DEL", "DUP", "INS", "INV")

for(a in 1:8){
  vcf <- readVcf(paste0("./Manta_20Mb_Masked_Results/", vcf.files[a]))
  #print(dim(vcf));
  GRDataFrame[[a]] <- as.data.frame(breakpointRanges(vcf))
  GRDataManta.tab <- cbind(GRDataManta.tab, as.data.frame(table(GRDataFrame[[a]]$svtype))[,2])
}

GRDataManta.tab <- GRDataManta.tab[,-1]
colnames(GRDataManta.tab) <- vcf.names
GRDprop.mat <- t(GRDataManta.tab)
GRDprop.mat1 <- GRDprop.mat/rowSums(GRDprop.mat)*100

write.csv(GRDprop.mat1, file = paste0("./Manta_20Mb_Masked_Results/CTR_edw23_0001-SVs_Manta_annoSummary_proportion_12_05_2019.csv"))
write.csv(GRDprop.mat, file = paste0("./Manta_20Mb_Masked_Results/CTR_edw23_0001-SVs_Manta_annoSummary_numbers_12_05_2019.csv"))

message("+------Produce Barplot for Manta proportion and VariantAnnotation proportion of SVs across samples----------+")

## format the data with sampleID, annotation, proportion for ggplot, and also order the proportion from greatest-smallest.
## Manta data reformatted
GRDplotdat <- as.data.frame(GRDprop.mat1)
GRDplotdat$ID <- rownames(GRDplotdat)
GRDplotdat.melt <- melt(GRDplotdat, id.var = 'ID')

GRDplotdat.melt <- within(GRDplotdat.melt, variable <- factor(variable, 
                                                              c('BND', 'INV', 'DUP', 'INS', 'DEL'), 
                                                              ordered = TRUE))
## Manta after VariantAnnotation reformatted
VRDplotdat <- as.data.frame(Freq.tableProp)
VRDplotdat$ID <- rownames(Freq.tableProp)
VRDplotdat.melt <- melt(VRDplotdat, id.var = 'ID')

VRDplotdat.melt <- within(VRDplotdat.melt, variable <- factor(variable, 
                                                              c('spliceSite', 'fiveUTR', 'threeUTR', 'coding', 'promoter', 'intergenic', 'intron'), 
                                                              ordered = TRUE))


message("+------- Split into control and mtrrgt for both SNPs and SVs -------------------+")

## Total SVs and SNPs annotation reformatted
SVs_avg <- colSums(Freq.tableProp)/8
SNPs_avg <- c(0.03, 58.296, 0.151, 0.279, 0.588, 40.619, 0.036)
names(SNPs_avg) <- names(SVs_avg)

SVs_SNPs.dat <- as.data.frame(rbind(SVs_avg, SNPs_avg))
rownames(SVs_SNPs.dat) <- c("SV", "SNP")
SVs_SNPs.dat$ID <- rownames(SVs_SNPs.dat) 
SVs_SNPs.melt <- melt(SVs_SNPs.dat, id.var = 'ID')
SVs_SNPs.melt <- within(SVs_SNPs.melt, variable <- factor(variable, 
                                                          c('spliceSite', 'fiveUTR', 'threeUTR', 'coding', 'promoter', 'intergenic', 'intron'), 
                                                          ordered = TRUE))

## SNPs annotation reformatted
SNP.samples <- matrix(0, nrow=8, ncol=7)
SNP.samples[1,] <- c(0.644,41.113,57.73,0.029,0.029,0.278,0.176)
SNP.samples[2,] <- c(0.598,41.102,57.799,0.028,0.035,0.292,0.146)
SNP.samples[3,] <- c(0.433,40.718,58.366,0.028,0.035,0.259,0.161)
SNP.samples[4,] <- c(0.451,40.99,58.013,0.029,0.036,0.298,0.182)
SNP.samples[5,] <- c(0.633,40.749,58.079,0.029,0.043,0.288,0.18)
SNP.samples[6,] <- c(0.477,42.023,56.929,0.036,0.043,0.311,0.181)
SNP.samples[7,] <- c(0.455,41.35,57.696,0.028,0.035,0.276,0.17)
SNP.samples[8,] <- c(0.64, 41.449,57.404,0.027,0.04,0.28,0.16)
rownames(SNP.samples) <- vcf.names
colnames(SNP.samples) <- c("coding", "intergenic", "intron", "spliceSite", "promoter", "threeUTR", "fiveUTR")
SNP.samples <- as.data.frame(SNP.samples)
SNP.samples$ID <- vcf.names
SNP_samples.melt <- melt(SNP.samples, id.var= 'ID')
SNP_samples.melt <- within(SNP_samples.melt, variable <- factor(variable, 
                                                                c('spliceSite', 'fiveUTR', 'threeUTR', 'coding', 'promoter', 'intergenic', 'intron'), 
                                                                ordered = TRUE))

message("+--- ggplot bar for ctrl and mtrrgt for snps and svs              ----------+")


SVs_avg.mtrr <- colSums(Freq.tableProp[1:6,])/6
SVs_avg.ctrl <- colSums(Freq.tableProp[7:8,])/2
SNPs_avg.mtrr <- colSums(SNP.samples[1:6,1:7])/6
SNPs_avg.ctrl <- colSums(SNP.samples[7:8,1:7])/2

SVs_avg.dat <- as.data.frame(rbind(SVs_avg.ctrl, SVs_avg.mtrr))
colnames(SVs_avg.dat) <- c("SpliceSiteRegion", "Intron", "fiveUTR", "threeUTR", "Exon", "Intergenic", "Promoter")
rownames(SVs_avg.dat) <- c("C57BL/6J", "Mtrrgt/gt")

SNPs_avg.dat <- as.data.frame(rbind(SNPs_avg.ctrl, SNPs_avg.mtrr))

colnames(SNPs_avg.dat) <- c("Exon", "Intergenic", "Intron", "SpliceSiteRegion",  "Promoter", "threeUTR", "fiveUTR")
rownames(SNPs_avg.dat) <- c("C57BL/6J", "Mtrrgt/gt")

SVs_avg.dat$ID <- rownames(SVs_avg.dat) 
SVs_avg.dat <- SVs_avg.dat[,c(2,6,7,5,4,3,1,8)]
SVs_avg.melt <- melt(SVs_avg.dat, id.var = 'ID')
SVs_avg.melt <- within(SVs_avg.melt, variable <- factor(variable, 
                                                        c("SpliceSiteRegion","fiveUTR","threeUTR", "Exon","Promoter", "Intergenic","Intron"), 
                                                        ordered = TRUE))
SNPs_avg.dat$ID <- rownames(SNPs_avg.dat) 
SNPs_avg.dat <- SNPs_avg.dat[,c(3,2,5,1,6,7,4,8)]
SNPs_avg.melt <- melt(SNPs_avg.dat, id.var = 'ID')
SNPs_avg.melt <- within(SNPs_avg.melt, variable <- factor(variable, 
                                                          c("SpliceSiteRegion","fiveUTR","threeUTR", "Exon","Promoter", "Intergenic","Intron"),
                                                          ordered = TRUE))

message("+-----                        Fig S3 (C, D)          ----------------------------+")

pdf("./Manta_20Mb_Masked_Results/CTR_edw23_0001-SVs_MantaAnno_RVAnno_SNPs_Proportion_CtrlvsMtrr_S3_C.pdf", height = 5, width = 12)

b1 <- ggplot(SVs_avg.melt, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  ylab("Annotation Proportion") +
  coord_flip() +
  geom_col() +
  scale_fill_manual(values = c(Intron="blue", Intergenic="black", Exon="purple",
                               threeUTR="pink", fiveUTR="darkgreen", Promoter="orange",
                               SpliceSiteRegion="brown")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("SVs")

b2 <- ggplot(SNPs_avg.melt, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  ylab("Annotation Proportion") +
  coord_flip() +
  geom_col() +
  scale_fill_manual(values = c(Intron="blue", Intergenic="black", Exon="purple",
                               threeUTR="pink", fiveUTR="darkgreen", Promoter="orange",
                               SpliceSiteRegion="brown")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("SNPs")

plot_grid(b1,b2, ncol = 1, nrow=2)

dev.off()

message("+ -----------------------------------------------------------------------------------------------+")
message("+----- Fig 1 (A)   129 SNPs/Indel vcf frequency across Mtrr region, 1kb adjacent window      ----+")
message("+ -----------------------------------------------------------------------------------------------+")

## calling two vcf file before filter and after filter
Re13.vcf <- read.table("recalibrated_variants_chr13.vcf", header = F)[,-c(8:9)] # 126404
MF.vcf <- read.table("recalibrated_vcf_filter1_segDup_simpRL9_hete3_RefH8D14.recode.vcf", header = F)[,-c(8:9)] ## 9494
MF13.vcf <- subset(MF.vcf, MF.vcf[,1]=="13") # 1243


start.locs <- 68560780; end.locs <- 68582149

nfwindow <- 20000000 ## flanking region size around Mtrr gene
swin <- 1000000  ## sliding window across the chrosmome



startw <- seq(start.locs-nfwindow, end.locs+nfwindow, swin)
endw <- c(startw[-1], end.locs+nfwindow)

CwinGR <- GRanges(seqnames = "13", IRanges(start=startw,end=endw),strand="+")
S129GR <- GRanges(seqnames = "13", IRanges(start=S129PFW.vcf[,2],end=S129PFW.vcf[,2]),strand="+")
SMFGR <- GRanges(seqnames = "13", IRanges(start=S129PFW.vcf[,2],end=S129PFW.vcf[,2]),strand="+")
Freq.count <- countOverlaps(CwinGR, S129GR,maxgap=-1L, minoverlap=0L,
                            type=c("any"), ignore.strand=FALSE)

Freq.countF <- Freq.count/sum(Freq.count)*100
xplot <- apply(cbind(startw, endw), 1, function(x) mean(x))
# add exon position into Mtrr region annotation
Mtrr_anno <- read.table(paste0(out.dir, "/Mtrr_Chr13_annot.txt"), header =F)
exons <- subset(Mtrr_anno, Mtrr_anno[,3]=="exon")


pdf(paste0(out.dir, "/CTR_edw23_0001-S129SNPs_FreqCounts_aw20Mb_sw1Mb_Mtrr.pdf"), width=7, height=5)
plot(xplot, Freq.countF, type="b", xlab="locs", ylab = "Frequnency of SNPs count/kb (%)")
abline(v=start.locs, lty=1, col = "red", lwd = 3)
abline(v=end.locs, lty=1, col ="red", lwd = 3)

dev.off()

message("+----------- Epigenetic memory of sperm DMRs Supp Figure 7, used Deeptools ----------+")
message("+-------Epigenetic memory of sperm DMRs Supp Figure 9, 10 used IGV  -----------------+")

##--------------------- FIN---------------------------------------------------------------------##

