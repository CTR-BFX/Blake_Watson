#!/usr/local/bin/Rscript
#
# Description
# Analysis R-Script for MEDIPS Template
#
# Author: Russell Hamilton (rsh46@cam.ac.uk )
# Date:   20170428 
#
#------------------------------------------------------------------------------


#source("http://bioconductor.org/biocLite.R")
#biocLite("MEDIPS")
library("MEDIPS")

biocLite("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Mmusculus.UCSC.mm10")

seqlengths(BSgenome.Mmusculus.UCSC.mm10)

# Set the current working directory where the bam files are
wd <- "/Users/rhamilto/Documents/CTR-Groups/Erica_Watson/GinaBlake_MeDiPs/"
setwd(wd)


# Some MeDIPs specific options

BSgenome   = "BSgenome.Mmusculus.UCSC.mm10"
uniq       = 1e-3 
extend     = 300
ws         = 100  #500
shift      = 0
chr.select = c("chr9")

chr.select = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
               "chr20", "chr21", "chr22")

# Load in the bam files "subdirectory" "bamfile"

# Group 1
Sample1 = "C57_35.chr.bam"
Sample2 = "sample2.bam"

# Group 2
Sample3 = "sample3.bam"
Sample4 = "sample4.bam"

BAM_File_List <- c( Sample1 )


#
# Saturation QC Plots 
#
for(i in 1:length(BAM_File_List)) 
{
  pdf(paste(wd,as.character(BAM_File_List[i]),"_saturationPlot.pdf",sep=""), width=7,height=7)
  sr=MEDIPS.saturation(file=BAM_File_List[i], BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
  MEDIPS.plotSaturation(sr)
  dev.off()
}


#
# Sequence pattern ("CG") coverage analysis to plot Pie charts and histograms
#
for(i in 1:length(BAM_File_List)) 
{
  cr = MEDIPS.seqCoverage(file = BAM_File_List[i], pattern = "CG", BSgenome = BSgenome, chr.select = chr.select, extend = extend, shift = shift, uniq = uniq)
  
  pdf(paste(wd,as.character(BAM_File_List[i]),"_coveragePie.pdf",sep=""), width=7,height=7)
  MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0, 5, 10, 20, 25, 50, 100, 250))
  dev.off()
  
  pdf(paste(wd,as.character(BAM_File_List[i]),"_coverageHist.pdf",sep=""), width=7,height=7)
  MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "hist", t = 100, main = "Sequence pattern coverage, histogram")
  dev.off()
}

#
# CpG Enrichment analysis
#
#er = MEDIPS.CpGenrich(file = ?????, BSgenome = BSgenome,chr.select = chr.select, extend = extend, shift = shift, uniq = uniq)


# MSET Group 1 
Group.1.MSET = MEDIPS.createSet(file = Sample1, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)
Group.1.MSET = c(Group.1.MSET, MEDIPS.createSet(file = Sample2, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws,chr.select = chr.select))

# MSET Group 1 
Group.2.MSET = MEDIPS.createSet(file = Sample3, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)
Group.2.MSET = c(Group.2.MSET, MEDIPS.createSet(file = Sample4, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws,chr.select = chr.select))


MSet_List = c(Group.1.MSET, Group.2.MSET)

# Calculate coupling vector for CpG density dependent normalization of MeDIP-seq data
CS.Group.1  = MEDIPS.couplingVector(pattern = "CG", refObj = Group.1.MSET)
CS.Group.2  = MEDIPS.couplingVector(pattern = "CG", refObj = Group.2.MSET)



#
# Export wiggle file for the coupling set
#
MEDIPS.exportWIG(Set = CS.Group.1, file = "CS.Group.1.wig", format = "pdensity", descr = "Coupling Plot CS.Group.1")


#
# Create wiggle files for MSets and coupling sets and Calculate and plot a calibration curve
#
for(i in 1:length(MSet_List)) 
{  
  MEDIPS.exportWIG(Set = MSet_List[i], file = paste(wd,as.character(MSet_List[i]),".wig",sep=""), format = "rpkm", descr = as.character(basename(MSet_List[i])))
  
  pdf(paste(wd,as.character(MSet_List[i]),"_callibrationPlot.pdf",sep=""),width=7,height=7)
  MEDIPS.plotCalibrationPlot(CSet = CS.OMNI.BL, main = as.character(basename(MSet_List[i])), MSet = MSet_List[i], plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
  dev.off()
}


#
# Differential methylation between two conditions
#
mr.edgeR = MEDIPS.meth(MSet1 = Group.1.MSET, MSet2 = Group.2.MSET, CSet = CS, p.adj = "bonferroni", diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1)

# Select significant windows
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1, adj = T, ratio = NULL, bg.counts = NULL, CNV = F)

dm1=mr.edgeR.s$edgeR.adj.p.value<0.1
dm2=mr.edgeR.s$edgeR.p.value<0.01
write.table(mr.edgeR.s[which(dm1),], paste(wd,"Group.1_Vs_Group.2_tab1_dmr.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

# Gain in Group.1 over Group.2
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] > 0), ]

mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.gain, distance=1)

dm1=mr.edgeR.s.gain.m$edgeR.adj.p.value<1
dm2=mr.edgeR.s.gain.m$edgeR.p.value<1
write.table(mr.edgeR.s.gain.m[which(dm1),], paste(wd,"Group.1_Vs_Group.2_tab1_dmr_gain.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)


# Loss in Group.1 over Group.2
mr.edgeR.s.loss = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] < 0), ]

mr.edgeR.s.loss.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.loss, distance=1)

dm1=mr.edgeR.s.loss.m$edgeR.adj.p.value<0.1
dm2=mr.edgeR.s.loss.m$edgeR.p.value<0.01
write.table(mr.edgeR.s.loss.m[which(dm1),], paste(wd,"Group.1_Vs_Group.2_tab1_dmr_loss.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)



#
# END OF SCRIPT
#
