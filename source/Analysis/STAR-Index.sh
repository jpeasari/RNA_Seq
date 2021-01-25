
#!/bin/bash 
mkdir -p ../../genome_index
cd ../../genome_index

STAR	 --runThreadN 12  \
       --runMode genomeGenerate \
       --genomeDir genome_index \
       --genomeFastaFiles ../../Data/raw/GCF_000146045.2_R64_genomic.fna \
       --sjdbGTFfile ../../Data/raw/GCF_000146045.2_R64_genomic.gff \
       --sjdbOverhang 49 \
       --genomeSAindexNbases 10  \
       --limitGenomeGenerateRAM 10000000000
