#!/bin/bash

echo " "
echo "Script to ensure pyTMB.py is working with SHIV01 WES vcfs."
echo " "

#TMB LOW sample

echo " "
echo "TMB Low sample :"
echo " "

python ../../bin/pyTMB.py \
    -i ./M1119-Tumor.tumor_constit.mutect.hg19_multianno.TUMOR_M1119-Tumor.Mutect.vcf.gz \
    --effGenomeSize 33280000\
    --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2 \
    --filterLowQual --filterNonCoding --filterSyn --filterPolym  \
    --polymDb 1k,gnomad  \
    --dbConfig ../../config/annovar_102015.yml \
    --varConfig ../../config/mutect2.yml

#TMB High sample

echo " "
echo "TMB High sample :"
echo " "

python ../../bin/pyTMB.py \
    -i ./MR293-Tumor.tumor_constit.mutect.hg19_multianno.TUMOR_MR293-Tumor.Mutect.vcf.gz \
    --effGenomeSize 33280000\
    --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2 \
    --filterLowQual --filterNonCoding --filterSyn --filterPolym  \
    --polymDb 1k,gnomad  \
    --dbConfig ../../config/annovar_102015.yml \
    --varConfig ../../config/mutect2.yml

echo " "
echo "Finished"
echo " "
