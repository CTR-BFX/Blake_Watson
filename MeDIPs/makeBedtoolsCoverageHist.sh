#!/bin/bash


GENOME_FILE="Mus_musculus.GRCm38.dna.chromosome.all.fa.chr.fai"


for i in *.srtd.bam;
do
  echo ${i}
  bedtools coverage -sorted -a Mus_musculus.GRCm38.100kbwindows.bed -b ${i} -g ${GENOME_FILE} > ${i/.bam/.100Kb.bed} 
  bedtools coverage -sorted -a Mus_musculus.GRCm38.25kbwindows.bed  -b ${i} -g ${GENOME_FILE} > ${i/.bam/.25Kb.bed}
  bedtools coverage -sorted -a Mus_musculus.GRCm38.5kbwindows.bed   -b ${i} -g ${GENOME_FILE} > ${i/.bam/.5Kb.bed}
  bedtools coverage -sorted -a Mus_musculus.GRCm38.500bpwindows.bed -b ${i} -g ${GENOME_FILE} > ${i/.bam/.500bp.bed}
done
