#!/bin/bash

echo " "
echo "Script to ensure pyTMB.py is working with SeqOIA WES vcfs."
echo " "


time python ../../bin/pyTMB.py \
    -i head_vcf_seqOIA.vcf.gz \
    --effGenomeSize 33958083 \
    --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 3 \
    --filterLowQual --filterNonCoding --filterPolym --filterSyn \
    --polymDb 1k,gnomad \
    --dbConfig ../../config/snpeff.yml \
    --varConfig ../../config/mutect2.yml
