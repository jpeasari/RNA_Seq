#!/bin/bash
mkdir -p ../../Data/feature_counts
cd ../../Data/feature_counts

featureCounts    -T 10 \
                 -t exon  \  
                 -g Parent  \ 
                 -a ../../Data/raw/R64_genomic_gff.gtf   \
                 -o ../../Data/feature_counts/feature_count.txt  \
                 ../star_align/*bam
