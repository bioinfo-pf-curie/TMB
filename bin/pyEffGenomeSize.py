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
    parser.add_argument("--featureTypes", help="List of Features to keep (column 3 from gtf/gff eg: exon, gene, Selenocysteine, start_codon, stop_codon, transcript, UTR, CDS", nargs='+', default = [])

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
    function using pybedtools nomenclature to filter transcript_types from gtf file with defined features.
    """
    #print(featuretype)
    if interval.attrs['transcript_type'] in featuretype:
        #print(interval.attrs['transcript_type'])
        return interval.attrs['transcript_type']
    return False

def filterFeatureGtf(interval, features):
    """
    function using pybedtools nomenclature to parse feature types such as exon, introns, CDS from a gtf file
    """
    if interval[2] in features:
        return interval[2]
    return False


    # intervalAnnot = interval.attrs.items()
    # #print(intervalAnnot)
    # print(interval.attrs['gene_type'], interval.attrs['transcript_type'])
    # for annot in interval.attrs['transcript_type']:
    #     print(annot)
    #     if annot in featuretype:
    #         return interval
    #return False

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
    with open(infile,'rt') as f:
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
    print("Loading data")
    myBed = pbt.BedTool(args.bed)
    myGtf = pbt.BedTool(args.gtf)

    coding = ["antisense", "IG_C_gene","IG_D_gene","IG_J_gene", "IG_V_gene","protein_coding","TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene","nonsense_mediated_decay","non_stop_decay"]
    nonCoding = ["3prime_overlapping_ncrna", "IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","lincRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","polymorphic_pseudogene","processed_transcript","processed_pseudogene","pseudogene","retained_intron","rRNA","sense_intronic","sense_overlapping","snoRNA","snRNA","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","unitary_pseudogene","unprocessed_pseudogene"]

    if args.filterCoding:
        filteredGtf = myGtf.filter(filterGtf, nonCoding).saveas("filtered_gtf.gtf")
        filteredGtf = pbt.BedTool(filteredGtf)
        #print(pbt.BedTool(filteredGtf).head())


    if args.filterNonCoding:
        filteredGtf = myGtf.filter(filterGtf, coding).saveas("filtered_gtf.gtf")
        filteredGtf = pbt.BedTool(filteredGtf)
        #print(pbt.BedTool(filteredGtf).head())

    if args.featureTypes:
        features = ["exon", "gene", "Selenocysteine", "start_codon", "stop_codon", "transcript", "UTR", "CDS"]
        for elem in args.featureTypes:
            if elem in features:
                featuretypes = args.featureTypes
                feature_filteredGtf = filteredGtf.filter(filterFeatureGtf, featuretypes).saveas("feature_filtered_gtf.gtf")
                filteredGtf = pbt.BedTool(feature_filteredGtf)
                print(pbt.BedTool(filteredGtf).head())
            else:
                print(type(args.featureTypes))
                print(type(features))
                sys.stderr.write("Error : wrong featureType. \nNot in [exon,gene,Selenocysteine,start_codon,stop_codon, transcript,UTR,CDS] list !?\n")
                sys.exit(-1)

        #print(pbt.BedTool(filteredGtf).head())




    ########################################
    ## mosdepth calculation
    if args.verbose:
        print("running mosdepth ...")
    sbp.run(["mosdepth", "-t", str(args.thread), "--by", str(args.bed), "-n", "--thresholds", str(args.minCoverage), "--mapq", str(args.minMapq), str(args.oprefix), str(args.bam)])

    # Unzip mosdepth bed file for intersect
    sbp.run(["gunzip", args.oprefix +".thresholds.bed.gz"])

    if args.verbose:
        print("running bedtools intersect ...")

    mosdepth_bed = pbt.BedTool(args.oprefix +".thresholds.bed")
    intersectBed = mosdepth_bed.intersect(filteredGtf).saveas('intersect_bed_gtf.bed')

    if args.verbose:
        print("head of intersect bed")
        print(intersectBed.head())

    effsize=getEffGenomeSizeFromMosdepth(args.oprefix +".thresholds.bed")
    effsize=getEffGenomeSizeFromMosdepth("intersect_bed_gtf.bed")

    ## clean modepth output files
    os.remove(args.oprefix + ".mosdepth.global.dist.txt")
    os.remove(args.oprefix + ".mosdepth.region.dist.txt")
    os.remove(args.oprefix + ".mosdepth.summary.txt")
    if not args.saveIntermediates:
        os.remove(args.oprefix + ".regions.bed.gz")
        os.remove(args.oprefix + ".thresholds.bed.gz")
    os.remove(args.oprefix + ".thresholds.bed.gz.csi")
    os.remove(args.oprefix + ".regions.bed.gz.csi")
