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

__version__ = '0.1.0dev'

"""
This script is designed to calculate an exact effective genome size according to ;
- The sample metrics (sequencing depth, mapping quality)
- The genomic regions considered for the variant calling (coding/non-coding regions)

python pyEffGenomeSize.py -bed ${BED} \
-bam ${BAM} --minCoverage 100 --minMapq 20 \
-gtf ${GTF} --filterNonCoding --filterCoding \
> result_folder
"""

import pybedtools as pbt
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
    parser.add_argument("--bed", help="BED file (.bed)", default=None)
    parser.add_argument("--gtf", help="GTF file for genome annotation (.gtf)", default=None)
    parser.add_argument("--bam", help="BAM file for mapping statistics (.bam)", default=None)

    #Filters
    parser.add_argument("--minCoverage", help="minimum coverage of the region",type=int, default=0)
    parser.add_argument("--minMapq", help="minimum coverage of the region",type=int, default=0)
    parser.add_argument("--filterNonCoding", help="Filter regions associated with non coding annotations",action="store_true")
    parser.add_argument("--filterCoding", help="Filter regions associated with coding annotations",action="store_true")

    # Others
    parser.add_argument("--saveIntermediates", help="Save mosdepth intermediate files",action="store_true")
    parser.add_argument("-t", "--thread", help="Number of threads for mosdepth run",type=int, default=1) 
    parser.add_argument("--oprefix", help="Suffix for filtered bed", type=str, default='pyeffg')
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s ("+__version__+")")

    args = parser.parse_args()
    return (args)

def filterGtf(interval, featuretype):
    """
    function using pybedtools nomenclature to parse gtf file with defined features.
    """
    intervalAnnot = interval.attrs.items()
    for annot in intervalAnnot:
        if annot[1] in featuretype:
            return interval
    return False

#def subsetFeaturetypes(featuretype):
#    """
#    Returns the filename containing only `featuretype` features.
#    """
#    return g.filter(filterGtf, featuretype).saveas().fn

def getEffGenomeSizeFromMosdepth(infile):
    """
    Calculate Effective Genome Size from a BED file
    """
    effgs = 0
    totgs = 0
    nline = 0
    with gzip.open(infile,'rt') as f:
        next(f)
        for line in f:
            bedtab = line.strip().split("\t")
            nline += 1
            try:
                chromosome, start, end, region, coverage = bedtab[:5]
            except ValueError:
                sys.stderr.write("Error : wrong input format in line", nline, ". Not a BED file !?")
                sys.exit(-1)

            intl = abs(int(end) - int(start))
            totgs += intl
            effgs += int(coverage)

    f.close()

    print("## Total region = {} (100%)".format(totgs))
    print("## Callable region = {} ({}%)".format(effgs, round(effgs/totgs*100, 3)))

    return effgs

if __name__ == "__main__":

    # Parse inputs
    args = argsParse()

    # Creating pybedtools objects
#    print("Loading data")
#    myBed = pbt.BedTool(args.bed)
#    myGtf = pbt.BedTool(args.gtf)

#    featuretypes = ['protein_coding', 'pseudogene']
#    coding = ["antisense", "IG_C_gene","IG_D_gene","IG_J_gene", "IG_V_gene","protein_coding","TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene"]
#    nonCoding = ["3prime_overlapping_ncrna", "IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","lincRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","polymorphic_pseudogene","processed_transcript","pseudogene","rRNA","sense_intronic","sense_overlapping","snoRNA","snRNA","TR_J_pseudogene","TR_V_pseudogene"]

#    if args.filterCoding:
#        filteredGtf = myGtf.filter(filterGtf, nonCoding).saveas("filtered_gtf.gtf")
#        print(pbt.BedTool(filteredGtf).head())

#    if args.filterNonCoding:
#        filteredGtf = myGtf.filter(filterGtf, coding).saveas("filtered_gtf.gtf")
#        filteredGtf = pbt.BedTool(filteredGtf)

#    print("Intersect")
#    intersectBed = myBed.intersect(filtered_gtf).saveas('intersect_bed_gtf.bed')
#    print(intersectBed.head())

    ########################################
    ## mosdepth calculation
    if args.verbose:
        println("running mosdepth ...")
    sbp.run(["mosdepth", "-t", str(args.thread), "--by", str(args.bed), "-n", "--thresholds", str(args.minCoverage), "--mapq", str(args.minMapq), str(args.oprefix), str(args.bam)])
    effsize=getEffGenomeSizeFromMosdepth(args.oprefix +".thresholds.bed.gz")
 
    ## clean modepth output files
    os.remove(args.oprefix + ".mosdepth.global.dist.txt")
    os.remove(args.oprefix + ".mosdepth.region.dist.txt")
    os.remove(args.oprefix + ".mosdepth.summary.txt")
    if not args.saveIntermediates:
        os.remove(args.oprefix + ".regions.bed.gz")
        os.remove(args.oprefix + ".thresholds.bed.gz")
    os.remove(args.oprefix + ".thresholds.bed.gz.csi")
    os.remove(args.oprefix + ".regions.bed.gz.csi")

