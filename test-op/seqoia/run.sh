#!/bin/bash

echo " "
echo "Script to ensure pyTMB.py is working with SeqOIA WES vcfs."
echo " "

SCRIPT_DIR=/bioinfo/users/tgutman/Documents/Tom/Pipelines/tmb;
INPUT_DIR=/data/tmp/tgutman/SeqOIA/TMB

time python ../../bin/pyTMB.py \
    -i $INPUT_DIR/MAP586_cancer_curie_final_no_reheader.vcf \
    --bed $INPUT_DIR/Merge_CORE-SPIKE_GRCh38.merged.bed \
    --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 3 \
    --filterLowQual --filterNonCoding --filterPolym --filterSyn \
    --polymDb 1k,gnomad \
    --dbConfig ../../config/snpeff.yml \
    --varConfig ../../config/mutect2.yml
