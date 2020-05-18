#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# Mask Mtrr region for DMRs
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/Blake_Watson
#
#
# Analysis Performed by Xiaohui Zhao
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

library(baySeq)


message("+---           Remove Mtrr 20Mb region                 ----------------+")

setwd("/storage/CTR-Projects/CTR_edw23/CTR_edw23_0001")
Chr.sizes.dat <- read.table("./CellReport_Jung_2017_Seq/mm10.chrom.sizes", header = F)
## Remove chrM, chrX and chrY
Chr.sizes <- Chr.sizes.dat[c(11:19,1:10),2]
Chrs <- Chr.sizes.dat[c(11:19,1:10),1]
dmr.names <- list.files(path="./PostGATK/DMRs_Data/", pattern = "*_allcol.bed")
dmr.files <- paste0("./PostGATK/DMRs_Data/", dmr.names)[1:3]
dmr.data <- lapply(dmr.files, function(x) read.table(x,header=F))
mtrr_mask20mb <- c(58060780,  80060780)
## only dmr.data +/gt and gt/gt has mtrr region dMRs

dmr.data[[1]] <- dmr.data[[1]][-which(dmr.data[[1]]=="chrUn_random"| dmr.data[[1]]=="chrY_random"), ]
dmr.data[[2]] <- dmr.data[[2]][-which(dmr.data[[2]]=="chrX"| dmr.data[[2]]=="chrY_random"), ]
dmr.data[[3]] <- dmr.data[[3]][-which(dmr.data[[3]]=="chrX"), ]

dmr.data[[2]] <- dmr.data[[2]][-which(dmr.data[[2]][,1]=="chr13"&dmr.data[[2]][,2] >mtrr_mask20mb[1]&dmr.data[[2]][,3]<mtrr_mask20mb[2]),]
dmr.data[[3]] <- dmr.data[[3]][-which(dmr.data[[3]][,1]=="chr13"&dmr.data[[3]][,2] >mtrr_mask20mb[1]&dmr.data[[3]][,3]<mtrr_mask20mb[2]),]

dmr.gr <- lapply(dmr.data, function(x) GRanges(seqnames=x[,1], IRanges(start=x[,2], end= x[,3])))
newdmr <- rbind(dmr.data[[1]], dmr.data[[2]], dmr.data[[3]])
newdmr.unique <- unique(newdmr)
newdmr <- newdmr[order(newdmr[,1], newdmr[,2]),]
write.table(newdmr, file = "./PostGATK/DMRs_Data/AllDMRs_maskMtrr20Mb_N459.bed", sep="\t", row.names=F, col.names=F, quote=F)

message("+---           FINISH                ----------------+")
