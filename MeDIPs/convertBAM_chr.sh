#!/bin/bash


for i in *.bam;
do
  echo "Extracting headers from ${i} and converting e.g. 1 to chr1"
  samtools view -H ${i} | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' > ${i/.bam/.header}
  echo "Applying the new chromosomes to reads in bam file"
  samtools reheader ${i/.bam/.header} ${i} > ${i/.bam/.chr.bam}
done
