#!/bin/bash
mkdir ../../star_align
cd ../../star_align

for sample in ../../Data/raw/*.fastq;
out = ${i}
do
  STAR --runMode alignReads \
        --runThreadN 20 \
        --genomeDir ../../genome_index \
        --outSAMtype BAM  SortedByCoordinate \
        --outSAMunmapped Within  \
        --outSAMattributes Standard \
        --outFileNamePrefix ../../star_align/out  \
        --readFilesIn ${sample} \
        --limitBAMsortRAM 10000000000
done
