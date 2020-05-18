#!/usr/local/bin/bash

#------------------------------------------------------------------------------
# TITLE
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

# Need a chromosome size file
# cat Mus_musculus.GRCm38.dna.chromosome.all.fa.fai | awk '{ print "chr"$1"\t"$2 } ' | sed 's/^chrMT/chrM/g' > /storage/Genomes/Mus_musculus/GRCm38/mm9.chr.new.sizes
# OR
# fetchChromSizes mm9 > /storage/Genomes/Mus_musculus/GRCm38/mm9.chrom.sizes

module load bedtools2/2.26.0
module load hgdownload.cse.ucsc.edu/20190415

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
       echo "Need 2 parameters! <bedgraph> <chrom info>"
     exit
     fi
     
     F=$1
     G=$2
     
     bedtools slop -i ${F} -g ${G} -b 0 | bedClip stdin ${G} ${F}.ori.clip
     
     LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s ${F}.ori.clip > ${F}.clip.tmp
     
     bedtools merge -i ${F}.clip.tmp -c 5 -d 0 -o mean > ${F}.clip
     
     LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n ${F}.clip > ${F}.sort.clip
     
     cut -f1,2,3,4 ${F}.sort.clip > ${F}.sel.sort.clip
     
     bedGraphToBigWig ${F}.sel.sort.clip ${G} ${F}.bw
     
     rm -f ${F}.ori.clip ${F}.clip.tmp ${F}.clip ${F}.sort.clip
 ## bash BedClip_MACS_overlap.sh TSC_H3K4me1_sorted.bed mm9.chrs.new.sizes