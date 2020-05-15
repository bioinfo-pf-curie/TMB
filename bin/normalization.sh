#!/bin/bash

#Normalization script with BCFTOOLS

while getopts "p:v:f:o:h" option; do
    case "${option}" in
        p) BCFTOOLS_DIR=${OPTARG};;
        v) VCF=${OPTARG};;
        f) FASTA=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        h)  echo "Usage:"
            echo "bash normalization.sh -p /path/to/bcftools/bin"
            echo "                      -v file.vcf"
            echo "                      -f reference.fa"
            exit
    esac
done

echo " "
echo "bcftools path :" $BCFTOOLS_DIR
echo "vcf file:" $VCF
echo "reference genome:" $FASTA
echo "output path": $OUTPUT
echo " "

#Count how many variants in the vcf and how many are multiallelic
NB_VAR_MULTI=`$BCFTOOLS_DIR/bcftools view -H -m 3 $VCF |wc -l`
NB_VAR=`$BCFTOOLS_DIR/bcftools view -H $VCF |wc -l`

echo " "
echo "Number of multiallelic variants"
echo $NB_VAR_MULTI
echo "Total number of variants"
echo $NB_VAR
echo " "

#Normalization step :

$BCFTOOLS_DIR/bcftools norm -f $FASTA -m- -o $OUTPUT/`basename $VCF .vcf`"_norm.vcf" $VCF

echo " "
echo "normalized vcf:" $OUTPUT/`basename $VCF .vcf`"_norm.vcf"
echo "finished"
echo " "
