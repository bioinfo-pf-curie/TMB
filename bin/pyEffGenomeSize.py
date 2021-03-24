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
--oprex /OUTPUT/PATH/suffix_file
"""

import argparse
import sys
import warnings
import os.path
import pandas as pd
import subprocess as sbp
import pybedtools as pbt

def argsParse():
    """
    Parse inputs
    """
    #Inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bed", help="BED file (.bed)", required=True, default=None)
    parser.add_argument("--gtf", help="GTF file for genome annotation (.gtf)", required=True, default=None)
    parser.add_argument("--bam", help="BAM file for mapping statistics (.bam)", default=None)

    #Filters
    parser.add_argument("--mosdepth", help="enable mosdepth processing, require .bam file", action="store_true")
    parser.add_argument("--minCoverage", help="minimum coverage of the region", type=int, default=0)
    parser.add_argument("--minMapq", help="minimum coverage of the region", type=int, default=0)
    parser.add_argument("--filterNonCoding", help="Filter regions associated with non coding annotations", action="store_true")
    parser.add_argument("--filterCoding", help="Filter regions associated with coding annotations", action="store_true")
    parser.add_argument("--featureTypes", help="List of Features (exon, gene, transcript, UTR, CDS) to keep (3rd column from gtf/gff). Required with --filterCoding argument ", nargs='+', default=[])

    # Others
    parser.add_argument("--saveIntermediates", help="Save mosdepth intermediate files", action="store_true")
    parser.add_argument("-t", "--thread", help="Number of threads for mosdepth run", type=int, default=1)
    parser.add_argument("--oprefix", help="Suffix for filtered bed", type=str, default='pyeffg')
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s ("+__version__+")")

    args = parser.parse_args()
    return args

def filterGtf(interval, featuretype):
    """
    function using pybedtools nomenclature to filter transcript_types from gtf file with defined features.
    """
    if interval.attrs['transcript_type'] in featuretype:
        return interval.attrs['transcript_type']
    return False

def filterFeatureGtf(interval, features):
    """
    function using pybedtools nomenclature to parse feature types such as exon, introns, CDS from a gtf file
    """
    if interval[2] in features:
        return interval[2]
    return False

def getEffGenomeSizeFromMosdepth(infile):
    """
    Calculate Effective Genome Size from a BED file
    """
    effgs = 0
    totgs = 0
    nline = 0
    with open(infile, 'rt') as f:
        if args.mosdepth:
            next(f)
        for line in f:
            bedtab = line.strip().split("\t")
            nline += 1
            if args.mosdepth:
                try:
                    chromosome, start, end, region, coverage = bedtab[:5]
                    effgs += int(coverage)
                except ValueError:
                    sys.stderr.write("Error : wrong input format in line", nline, ". Not a BED file !?")
                    sys.exit(-1)

            else:
                try:
                    chromosome, start, end = bedtab[:3]
                except ValueError:
                    sys.stderr.write("Error : wrong input format in line", nline, ". Not a BED file !?")
                    sys.exit(-1)

            intl = abs(int(end) - int(start))
            totgs += intl
    f.close()
    print("## File Name: {}".format(infile))
    print("## Total region = {}".format(totgs))
    if args.mosdepth:
        print("## Callable region = {} ({}%)\n".format(effgs, round(effgs/totgs*100, 3)))

if __name__ == "__main__":

    # Parse inputs
    args = argsParse()

    # Warnings:
    if (args.bam or args.minCoverage or args.minMapq) and not args.mosdepth:
        sys.stderr.write("Error: --mosdepth is required if --bam, --minCoverage or --minMapq is used.\n")
        sys.exit(-1)

    # Creating pybedtools objects
    if args.verbose:
        print("[RUNNING INFO]: Loading data\n")
    myBed = pbt.BedTool(args.bed)
    myGtf = pbt.BedTool(args.gtf)

    # Filtering step:
    # List of annotation from transcript_type field in gtf:
    coding = ["antisense", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "protein_coding", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "nonsense_mediated_decay", "non_stop_decay"]
    nonCoding = ["3prime_overlapping_ncrna", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "lincRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "polymorphic_pseudogene", "processed_transcript", "processed_pseudogene", "pseudogene", "retained_intron", "rRNA", "sense_intronic", "sense_overlapping", "snoRNA", "snRNA", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"]

    if args.verbose:
        print("[RUNNING INFO]: filtering GTF file\n")

    if args.filterNonCoding:
        filteredGtf = myGtf.filter(filterGtf, coding)
        filteredGtf = pbt.BedTool(filteredGtf)
        featuretypes = ["exon"]
        feature_filteredGtf = filteredGtf.filter(filterFeatureGtf, featuretypes).saveas("filtered_gtf.gtf")
        filteredGtf = pbt.BedTool(feature_filteredGtf)
        print(pbt.BedTool(filteredGtf).head())

    if args.filterCoding:
        if not args.featureTypes:
            sys.stderr.write("Error: --filterCoding requires --featureTypes.\n")
            sys.exit(-1)
        else:
            filteredGtf = myGtf.filter(filterGtf, nonCoding)
            filteredGtf = pbt.BedTool(filteredGtf)

            features = ["exon", "gene", "transcript", "UTR", "CDS"]
            check = set(args.featureTypes).issubset(features)

            if check is True:
                featuretypes = args.featureTypes
                feature_filteredGtf = filteredGtf.filter(filterFeatureGtf, featuretypes).saveas("filtered_gtf.gtf")
                filteredGtf = pbt.BedTool(feature_filteredGtf)
                print(pbt.BedTool(filteredGtf).head())

            else:
                sys.stderr.write("Error : wrong featureType. \nNot in [exon,gene,transcript,UTR,CDS] list !?\n")
                sys.exit(-1)

    if args.filterNonCoding or args.filterCoding:
        if args.verbose:
            print("[RUNNING INFO]: running bedtools intersect on gtf and bed...\n")

        x = pbt.BedTool()
        intersectBed = x.multi_intersect(i=[myBed.fn, filteredGtf.fn], header=True)
        print(pbt.BedTool(intersectBed).head())
        df = pd.read_table(intersectBed.fn, low_memory=False)

        #Obtain features present in filtered gtf and bed file
        df = df[df.num == 2]
        df.columns = df.iloc[0]
        df = df.iloc[:, [0, 1, 2]]
        df[1:].to_csv(args.oprefix +".intersect.bed", sep='\t', index=False)

    if args.mosdepth:
        # Mosdepth calculation
        if args.verbose:
            print("[RUNNING INFO]: running mosdepth ...\n")

        sbp.run(["mosdepth", "-t", str(args.thread), "--by", str(args.oprefix +".intersect.bed"), "-n", "--thresholds", str(args.minCoverage), "--mapq", str(args.minMapq), str(args.oprefix), str(args.bam)])
        # Unzip mosdepth bed file
        sbp.run(["gunzip", args.oprefix +".thresholds.bed.gz"])
        getEffGenomeSizeFromMosdepth(args.oprefix +".thresholds.bed")

        ## clean modepth output files
        os.remove(args.oprefix + ".mosdepth.global.dist.txt")
        os.remove(args.oprefix + ".mosdepth.region.dist.txt")
        os.remove(args.oprefix + ".mosdepth.summary.txt")
        if not args.saveIntermediates:
            os.remove(args.oprefix + ".regions.bed.gz")
            os.remove(args.oprefix + ".thresholds.bed")
        os.remove(args.oprefix + ".thresholds.bed.gz.csi")
        os.remove(args.oprefix + ".regions.bed.gz.csi")

    else:
        getEffGenomeSizeFromMosdepth(args.oprefix +".intersect.bed")

    #END
    print("\nScript Ended Successfully !\n")
