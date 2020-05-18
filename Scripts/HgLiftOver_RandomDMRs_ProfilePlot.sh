# !/usr/bin/bash
#------------------------------------------------------------------------------
# liftover for public data mm9tomm10, and deeptools for profiles plot        
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

module load bedtools2/2.26.0
module load hgdownload.cse.ucsc.edu/20190415

## liftover

basedir=/storage/CTR-Projects/CTR_edw23/CTR_edw23_0001
mm9dir=/storage/CTR-Projects/CTR_edw23/CTR_edw23_0001/CellReport_Jung_2017_Seq
dmrdir=/storage/CTR-Projects/CTR_edw23/CTR_edw23_0001/PostGATK/DMRs_Data

cd ${mm9dir}
for names in *_sorted.bed 
do

  bash BedClip_MACS_overlap.sh ${names} mm9.chrs.new.sizes
  liftOver ${names}.bw mm9ToMm10.chain ${names}.mm10.bw
  
done

## Generate random Mtrr regions across chromosomes.

cd ${dmrdir}
bedtools shuffle -i AllDMRs_maskMtrr20Mb_N459.bed \
                 -g 550 Mus_GRCm38_chrall.fa.fai \
                 -chrom -seed 27442958 -noOverlapping -excl Mtrr_mask20Mb.bed \
                 > AllDMRs_maskMtrr20Mb_N459_random.bed

## compute the matrix and plot the profile 

cd ${basedir}


computeMatrix scale-regions -S ./ATAC-Seq/GSM2229960_Epi_Pair1_S1.bw \
./ATAC-Seq/GSM2229962_Exe_Pair1_S2.bw \
-R  ./PostGATK/DMRs_Data/New_Merge_C57_Vs_++_0.01_xy_6518_tab1_dmr_allcol_mm9.bed \
./PostGATK/DMRs_Data/Random_C57_Vs_++_0.01_xy_6518_tab1_dmr_allcol_mm9.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 3000 \
--afterRegionStartLength 5000 \
--skipZeros -o ./ATAC-Seq/Mtrr_C57_++_0.01_matrix_R3kb.mat.gz

computeMatrix scale-regions -S ./ATAC-Seq/GSM2229960_Epi_Pair1_S1.bw \
./ATAC-Seq/GSM2229962_Exe_Pair1_S2.bw \
-R  ./PostGATK/DMRs_Data/New_Merge_C57_Vs_+gt_0.01_tab1_dmr_xy_29418_allcol_mm9.bed \
./PostGATK/DMRs_Data/Random_C57_Vs_+gt_0.01_tab1_dmr_xy_29418_allcol_mm9.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 3000 \
--afterRegionStartLength 5000 \
--skipZeros -o ./ATAC-Seq/Mtrr_C57_+gt_0.01_matrix_R3kb.mat.gz

computeMatrix scale-regions -S ./ATAC-Seq/GSM2229960_Epi_Pair1_S1.bw \
./ATAC-Seq/GSM2229962_Exe_Pair1_S2.bw \
-R  ./PostGATK/DMRs_Data/New_Merge_C57_Vs_gtgt_0.01_xy_allcol_mm9.bed \
./PostGATK/DMRs_Data/Random_C57_Vs_gtgt_0.01_xy_allcol_mm9.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 3000 \
--afterRegionStartLength 5000 \
--skipZeros -o ./ATAC-Seq/Mtrr_C57_gtgt_0.01_matrix_R3kb.mat.gz


computeMatrix scale-regions -S ./ATAC-Seq/GSM2229960_Epi_Pair1_S1.bw \
./ATAC-Seq/GSM2229962_Exe_Pair1_S2.bw \
-R  ./PostGATK/DMRs_Data/New_Merge_C57_Vs_++_0.01_xy_6518_tab1_dmr_allcol_mm9.bed \
./PostGATK/DMRs_Data/New_Merge_C57_Vs_+gt_0.01_tab1_dmr_xy_29418_allcol_mm9.bed \
./PostGATK/DMRs_Data/New_Merge_C57_Vs_gtgt_0.01_xy_allcol_mm9.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 3000 \
--afterRegionStartLength 5000 \
--skipZeros -o ./ATAC-Seq/Mtrr_C57_++_+gt_gtgt_0.01_matrix_R3kb.mat.gz

computeMatrix scale-regions -S ./ATAC-Seq/GSM2229960_Epi_Pair1_S1.bw \
./ATAC-Seq/GSM2229962_Exe_Pair1_S2.bw \
-R  ./PostGATK/DMRs_Data/New_Merge_allDMRs_mm9.bed \
./PostGATK/DMRs_Data/Random_allDMRs_mm9.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 3000 \
--afterRegionStartLength 5000 \
--skipZeros -o ./ATAC-Seq/allDMRs_random_matrix_R3kb.mat.gz





plotHeatmap -m ./ATAC-Seq/Mtrr_C57_++_0.01_matrix_R3kb.mat.gz \
--colorMap Purples Greens \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs distance" \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--plotTitle "ATAC-Seq" \
--legendLocation none \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "Mtrr +/+" "Random" \
-out ./ATAC-Seq/Mtrr_C57_++_0.01_random_R3kb.png 

plotHeatmap -m ./ATAC-Seq/Mtrr_C57_+gt_0.01_matrix_R3kb.mat.gz \
--colorMap Purples Greens \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs distance" \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--plotTitle "ATAC-Seq" \
--legendLocation none \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "Mtrr +/gt" "Random" \
-out ./ATAC-Seq/Mtrr_C57_+gt_0.01_random_R3kb.png 

plotHeatmap -m ./ATAC-Seq/Mtrr_C57_gtgt_0.01_matrix_R3kb.mat.gz \
--colorMap Purples Greens \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs distance" \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--plotTitle "ATAC-Seq" \
--legendLocation none \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "Mtrr gt/gt" "Random" \
-out ./ATAC-Seq/Mtrr_C57_gtgt_0.01_random_R3kb.png 


plotHeatmap -m ./ATAC-Seq/Mtrr_C57_++_+gt_gtgt_0.01_matrix_R3kb.mat.gz \
--colorMap Purples Greens \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs distance" \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--plotTitle "ATAC-Seq" \
--legendLocation none \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "Mtrr +/+" "Mtrr +/gt" "Mtrr gt/gt" \
-out ./ATAC-Seq/Mtrr_C57_++_+gt_gtgt_0.01_R3kb.png 


## Merge all DMRs

plotHeatmap -m ./ATAC-Seq/allDMRs_random_matrix_R3kb.mat.gz \
--colorMap Purples Greens \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs distance" \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--plotTitle "ATAC-Seq" \
--legendLocation none \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "Mtrr all" "Random" \
-out ./ATAC-Seq/allDMRs_random_R3kb.png 


plotProfile -m matrix.mat.gz \
-out ExampleProfile2.png \
--colors red yellow \
--numPlotsPerRow 2 \
--yAxisLabel "DMRs" \
--startLabel Start \
--endLabel End \
--samplesLabel "Epiblast" "Extraembryonic" \
--regionsLabel "label 1" "label 2" \
--plotTitle "ATAC-Seq"


## histone modification from sperm data


computeMatrix reference-point -S ./CellReport_Jung_2017_Seq/GSM2088391_SPERM_H3K4me3.wig.bedgraph.bw \
-R  ./PostGATK/DMRs_Data/AllDMRs_maskMtrr20Mb.bed \
--referencePoint TSS \
--beforeRegionStartLength 1500 \
--afterRegionStartLength 1500 \
--binSize 10 \
--sortRegions descend \
--sortUsing mean \
--skipZeros -o ./CellReport_Jung_2017_Seq/AllDMRs_maskMtrr20Mb_H3K4me3.mat.gz

plotHeatmap -m ./CellReport_Jung_2017_Seq/AllDMRs_maskMtrr20Mb_H3K4me3.mat.gz \
--boxAroundHeatmaps yes \
--sortUsing mean \
--averageTypeSummaryPlot mean \
--heatmapWidth 10 \
--heatmapHeight 18 \
--xAxisLabel "DMRs start" \
--yAxisLabel "DMRs" \
--startLabel Start \
--plotTitle "H3K4me3-Seq" \
--legendLocation none \
-out ./CellReport_Jung_2017_Seq/AllDMRs_maskMtrr20Mb_H3K4me3.png 

#------------------FINIISH-------------------------------------------------

