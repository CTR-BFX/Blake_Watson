#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# WGS filter locs analysis, after first step filter1 by vcftools      
# Data from UCSC and aligned vcf files
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


message("+---------- Basic settings and libraries loadings --------------------+")

setwd("/Users/xz289/Documents/CTR_edw23_0001/UCSC_Data/")
source("/Users/xz289/Documents/CTR_edw23_0001/Scripts/Utility_functions.R")
library("segmentSeq") 
library("dplyr")
library("baySeq")

message("+---------- seqkit to get homo > 8 and dinu > 14 --------------------+")

system("seqkit locate -P -f Homo_Dinu_motif_all.fa Mus_musculus.GRCm38.dna.chromosome.all.fa > Filter_homoL8_dinuL14.bed")

message("+---- Useful functions for the filtering ------------------------------+")

annotation_fn <- function(annotations, chr.col, start.col, end.col, anno.col, outfile){
  anno.gtf <- read.table(annotations, sep="\t", quote="")
  anno.gtf.gene <- anno.gtf[anno.gtf[,anno.col]=="gene",]
  newanno <- anno.gtf[, c(chr.col, start.col, end.col, anno.col)]
  newanno.title <- c("#CHROM  START END Annotation")
  write.table(newanno.title, file = outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(newanno, file = outfile, append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#-----------------------------------------------------------------------------------------#
#     filtering based on UCSC Table Browser output                                        #
#-----------------------------------------------------------------------------------------#

filter_fn_bcf <- function(input.dat, inchr.col, instart.col, inwidth, bcfwin,
                          filter.dat, filchr.col, filstart.col, filend.col, 
                          output.dat){
  
  input <- read.table(input.dat, header = F)
  simpleR <- read.table(filter.dat, header = F)
  simpleRU <- distinct(simpleR, simpleR[,filchr.col], simpleR[,filstart.col], 
                       simpleR[,filend.col])
  
  # note the simpleGR start and end with +1 as in the UCSC locs start from 0.
  # bcfwin for if UCSC filter table is 1, otherwise is 0.
  inputGR <- GRanges(seqnames = paste0("chr", input[,inchr.col]), 
                     IRanges(start = input[,instart.col], width=1))
  
  simpleGR <- GRanges(seqnames = simpleRU[,filchr.col], 
                      IRanges(start = simpleRU[,filstart.col]+bcfwin, end = simpleRU[,filend.col]+bcfwin))
  
  simpleOver <- findOverlaps(inputGR, simpleGR)
  
  filterLoc <- as.data.frame(inputGR[queryHits(simpleOver)])[,1:2]
  filterLocU <- distinct(filterLoc, seqnames, start)
  filterLocU$seqnames  <- substr(filterLocU$seqnames, 4, 5)
  write.table(filterLocU, file = output.dat, append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  gc()
}

#-----------------------------------------------------------------------------------------------#
#  filter from seqkit output for homo and dinu                                                  #
#  we used Oey(2015) paper                                                                      #
#-----------------------------------------------------------------------------------------------#

filter_fn_seqkit <- function(input.dat, inchr.col, instart.col, inwidth, bcfwin,
                             filter.dat, filchr.col, filstart.col, filend.col, 
                             output.dat){
  
  input <- read.table(input.dat, header = F)
  simpleR <- read.table(filter.dat, header = T)
  simpleRU <- distinct(simpleR, simpleR[,filchr.col], simpleR[,filstart.col], 
                       simpleR[,filend.col])
  
  # note the simpleGR start and end with +1 as in the UCSC locs start from 0.
  # bcfwin for if UCSC filter table is 1, otherwise is 0.
  inputGR <- GRanges(seqnames = input[,inchr.col], 
                     IRanges(start = input[,instart.col], width=1))
  simpleGR <- GRanges(seqnames = simpleRU[,filchr.col], 
                      IRanges(start = simpleRU[,filstart.col]+bcfwin, end = simpleRU[,filend.col]+bcfwin))
  
  simpleOver <- findOverlaps(inputGR, simpleGR)
  
  filterLoc <- as.data.frame(inputGR[queryHits(simpleOver)])[,1:2]
  filterLocU <- distinct(filterLoc, seqnames, start)
  write.table(filterLocU, file = output.dat, append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  gc()
}

#-----------------------------------------------------------------------------------------------#
#  count heterzygosity within a given window                                                    #
#  we used Oey(2015) paper window within 10kb with more than 3 variants hete                    #
#-----------------------------------------------------------------------------------------------#

Heter_Wcount_fn <- function(input.dat, chr.col, locs.col, windowS, widthR, slideR, chrL.dat, output.dat){
  input <- read.table(input.dat, header = F)  # vcf file 
  chrL <- read.table(chrL.dat, header = F)  # fai file with each chr length
  chrL <- chrL[chrL[,chr.col]!="MT",]
  chrL.names <- chrL[,chr.col]
  finalLocs <- NULL
  for(chr in chrL.names){
    print(chr);
    chr.sub <- input[input[,chr.col]==chr,]
    chr.diff <- c(chr.sub[1,locs.col],diff(chr.sub[,locs.col]))
    chrdiff.ind <- ifelse(chr.diff > windowS, 0, 1)
    results <- data.frame(index = numeric(), win.sum = numeric())
    indexi <- 1;  # indexi the first value of the window for sliding
    indexj <- 1; # indexj the row of the results to be added next
    while(indexi < length(chrdiff.ind)) {
      #This is for window define calculate the sum.
      win.sum <- sum(chrdiff.ind[indexi:(indexi+widthR)], na.rm = TRUE)
      #Insert the results
      results[indexj, ] <- c(indexi, win.sum)
      #Increment the indices for the next pass
      indexi <- indexi + slideR
      indexj <- indexj + 1
    }
    hetewin.ind <- which(results$win.sum > 3)
    if(length(hetewin.ind)!=0){
      subhet.res <- results[hetewin.ind,]
      sublocs1 <- unique(unlist(apply(subhet.res, 1, function(x) list(c(x[1]:(x[1]+x[2]-1))))))
      chrlocs <- cbind(chr, chr.sub[sublocs1,2])
      colnames(chrlocs) <- c("chr", "locs")
      finalLocs <- rbind(finalLocs, chrlocs)
    }
    finalLocs
  }
  write.table(finalLocs, file = output.dat, append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  gc()
}

message("+----- Calling the UCSC filter files with simply repeats and segment duplications----+")

filter.dat <- "Filter_homoL8_dinuL14_simp.bed"
filter.dat1 <- "SimpleRepL9_Filter.txt"
filter.dat2 <- "SegmentDup_Filter.txt"
chrL.dat <- "Mus_musculus.GRCm38.dna.chromosome.all.fa.fai"

message("+----- calling the vcf filter file and remove the above UCSC filter files step by step----+")

MFinput <- "/Users/xz289/Documents/CTR_edw23_0001/VCF_Data/recalibrated_vcf_filter1.vcf"
MFoutputs <- paste0("recalibrated", c("_seqkit_HomDinuRef", "_simpRL9", "_segmeDup", "_heter3"), ".txt")

MF1 <- filter_fn_seqkit(MFinput, 1, 2, 1, 0, filter.dat, 1, 2, 3, MFoutputs[1])
MF2 <- filter_fn_bcf(MFinput, 1, 2, 1, 1, filter.dat1, 1, 2, 3, MFoutputs[2])
MF3 <- filter_fn_bcf(MFinput, 1, 2, 1, 1, filter.dat2, 1, 2, 3, MFoutputs[3])
MF4 <- Heter_Wcount_fn(MFinput, 1, 2, 10000, 4, 1, chrL.dat, MFoutputs[4])

system(paste("cat ", "recalibrated_segmeDup.txt ", "recalibrated_seqkit_HomDinuRef.txt ",  "recalibrated_simpRL9.txt |",
             "sort -k1 -k2 | uniq > ",  "recalibrated_exclude_locs.txt", sep=""))

message("+-----                FINISH                 ----+")

