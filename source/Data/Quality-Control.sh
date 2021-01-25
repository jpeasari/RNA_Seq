#!/bin/bash
mkdir -p ../../Reports/QC
cd ../../Reports/QC


for file in raw/*;
do
  fastqc -o = .$file
done

multiqc .
