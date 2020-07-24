#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of pyTMB software.
#
#  Copyright (c) 2020 - Institut Curie
#
#  File author(s):
#      Tom Gutman <tom.gutman@curie.fr>,
#      Nicolas Servant <nicolas.servant@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################

__version__ = '0.1.0'

"""
This script is designed to provide an ajusted bed file with it's corresponding size in order to be used for TMB calculation.

2 filters will be applied to the bed file provided by the user.
First a coverage filter will be applied. The treshold coverage will be defined by the user.
Then a filter based on the types of variants is also proposed. The user can specify what types of annotation is present in the region of the filtered bed.

python pyEffGenomeSize.py -bed ${BED} -bam ${BAM} -gtf ${GTF} \
-- filterNonCoding \
-- filterCoding \
-- minCoverage 100 > result_folder
"""

import pybedtools as pbt
from cyvcf2 import VCF
import argparse
import sys
import warnings
import re
import yaml
import os.path
import numpy as np
from datetime import date
import multiprocessing
import subprocess as sbp
import gzip

def argsParse():
    """
    Parse inputs
    """
    #Inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bed", help="Bed file to filter (.bed)", default=None)
    parser.add_argument("--gtf", help="gtf file for annotation (.gtf)", default=None)
    parser.add_argument("--bam", help="bam file to compute coverage (.bam)", default=None)

    #Filters
    parser.add_argument("--filterNonCoding", help="Filter regions associated with non coding annotations",action="store_true")
    parser.add_argument("--filterCoding", help="Filter regions associated with coding annotations",action="store_true")

    # Others
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s ("+__version__+")")
    parser.add_argument("-t", "--thread", help="Number of threads",type=int, default=1)
    parser.add_argument("--minCoverage", help="minimum coverage of the region",type=int, default=100)
    parser.add_argument("--outputName", help="Suffix for filtered bed",type=str)

    args = parser.parse_args()
    return (args)

def filter_gtf(interval, featuretype):
    """
    function using pybedtools nomenclature to parse gtf file with defined features.
    """
    interval_annot = interval.attrs.items()
    for annot in interval_annot:
        if annot[1] in featuretype:
            return interval
    return False

def subset_featuretypes(featuretype):
    """
    Returns the filename containing only `featuretype` features.
    """
    return g.filter(filter_gtf, featuretype).saveas().fn

def getEffGenomeSizeFromBed(infile, verbose=False):
    """
    Calculate Effective Genome Size from a BED file
    """
    if verbose:
        print("## Loading BED file '", infile, "'...")

    bedhandle = open(infile)
    effgs = 0
    nline = 0
    for line in bedhandle:
        bedtab = line.strip().split("\t")
        nline += 1
        try:
            chromosome, start, end = bedtab[:3]
        except ValueError:
            sys.stderr.write("Error : wrong input format in line", nline, ". Not a BED file !?")
            sys.exit(-1)

        intl = abs(int(end) - int(start))
        effgs += intl

    bedhandle.close()
    return effgs

if __name__ == "__main__":

    #Parse inputs
    args = argsParse()

    #Creating pybedtools objects
    print("Creating pybedtools objects")
    my_bed = pbt.BedTool(args.bed)
    my_gtf = pbt.BedTool(args.gtf)

    print("Filtering gtf")
    #g = pbt.BedTool(my_gtf)
    #pool = multiprocessing.Pool(args.thread)

    featuretypes = ['protein_coding', 'pseudogene']
    coding = ["antisense", "IG_C_gene","IG_D_gene","IG_J_gene", "IG_V_gene","protein_coding","TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene"]
    non_coding = ["3prime_overlapping_ncrna", "IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","lincRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","polymorphic_pseudogene","processed_transcript","pseudogene","rRNA","sense_intronic","sense_overlapping","snoRNA","snRNA""TR_J_pseudogene","TR_V_pseudogene"]

    if args.filterCoding:

        filtered_gtf = my_gtf.filter(filter_gtf, non_coding).saveas("filtered_gtf.gtf")
        print(pbt.BedTool(filtered_gtf).head())

    if args.filterNonCoding:

        filtered_gtf = my_gtf.filter(filter_gtf, coding).saveas("filtered_gtf.gtf")
        filtered_gtf = pbt.BedTool(filtered_gtf)

    print("Intersect")
    intersect_bed = my_bed.intersect(filtered_gtf).saveas('intersect_bed_gtf.bed')
    print(intersect_bed.head())

    threads, coverage, intersect_bed = str(args.thread), str(args.minCoverage), str(intersect_bed.fn)
    print(threads, coverage)

    print("Running mosdepth")
    #Bash command to calculate depth of coverage:
    sbp.run(["mosdepth", "-t", threads, "-b", intersect_bed, "-n", "-T", coverage, args.outputName, args.bam])

    #Genome size
    sbp.run(["gunzip", args.outputName +".regions.bed.gz"])

    effGS_intersect = getEffGenomeSizeFromBed(args.outputName +".regions.bed")
    #effGS_gtf = getEffGenomeSizeFromBed(args.gtf)
    effGS_bed = getEffGenomeSizeFromBed(args.bed)
    print(effGS_intersect, effGS_bed)

    #with gzip.open(args.outputName +".regions.bed.gz", 'rb') as bed:
        # print(bed)
        # unzip_bed = gzip.decompress(bed.read())
        # print(unzip_bed, bed)
        #effGS = getEffGenomeSizeFromBed(unzip_bed)
        #print(effGS)

#print number of lines with mean coverage >100
#zcat new_test.regions.bed.gz |awk  '{if($5>100)print $0}' |wc -l
#tail -n +6 /data/annotations/pipelines/Human/hg19/gtf/gencode.v19.annotation.sorted.gtf |cut -f9 |cut -d ";" -f3 |sort |uniq -c

#Running the tool:
#python pyEffGenomeSize.py --bed SureSelect_Clinical_Research_Exome_Regions_v2.bed --gtf head_gencode.gtf  --bam test_data/Preprocess/D251E01_uniq.nodup.onTarget.q20.recal.bam --thread 4 --minCoverage 150 --outputName test --filterNonCoding

#Problème à résoudre
- le filtre sur le coverage n'est pas fonctionnel,
- vérifier que tous les arguments sont bien utilisés (--verbose..)
- mieux gérer les filtres filterNonCoding et filterCoding, car pour le moment c'est soit l'un soit l'autre...
- bien tester le script pour s'assurer de son fonctionnement
- mettre des messages d'erreur si argument incorrect (try, except, sys.exit()...)
- mieux gérer les fichiers de sortie
- mettre à jour le yaml
