#!/bin/bash

echo "Running test on WES vcf files ..."

for i in $(ls *vcf.gz); do
    echo "$i ... "
    python ../../bin/pyTMB.py \
	   -i $i \
	   --effGenomeSize 33280000 \
	   --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2 \
	   --filterLowQual --filterNonCoding --filterSyn --filterPolym  \
	   --polymDb 1k,gnomad  \
	   --dbConfig ../../config/annovar_102015.yml \
	   --varConfig ../../config/mutect2.yml >> test.log
    echo "ok"
done

echo "done"
