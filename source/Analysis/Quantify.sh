#!/bin/bash
mkdir -p ../../Data/feature_counts
cd ../../Data/feature_counts

featureCounts    -T 10 \
                 -t exon  \  
                 -g Parent  \ 
                 -a ../../Data/raw/GCF_000146045.2_R64_genomic.gtf   \
                 -o ../../Data/feature_counts/feature_count.txt  \
                 ../../Data/raw/star_align/*bam
