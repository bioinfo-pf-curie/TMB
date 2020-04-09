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

import cyvcf2 
#import VCF
import argparse
import sys
import warnings
import re
import yaml
from datetime import date


def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")
    parser.add_argument("-r", "--rec", help="Recurence file")
    parser.add_argument("-o", "--out", help="Output file", default="output.vcf")
    args = parser.parse_args()
    return (args)


def loadRec(infile):

    handle = open(infile)
    nrec = 0
    nline = 0
    rec = {}
    pattern=re.compile("Recurrent")
    for line in handle:
        nline+=1
        bedtab = line.strip().split("\t")
        if re.search("Recurrent", bedtab[108]):
            k=str(bedtab[2])+":"+str(bedtab[3])+"-"+str(bedtab[4])
            rec[k] = bedtab[107]
            nrec +=1

    handle.close()
    return(rec)


if __name__ == "__main__":

    args = argsParse()

    ## Loading Data
    vcf = cyvcf2.VCF(args.vcf)
    rec = loadRec(args.rec)

    ## rec file
    vcf.add_info_to_header({'ID': 'RUNREC', 'Description': 'Run recurrence',                                                                                                                             
                                'Type':'Character', 'Number': '1'})                                                                                                                                      
    w = cyvcf2.Writer(args.out, vcf)

    for variant in vcf:
        k = str(variant.CHROM) + ":" + str(variant.start+1) + "-" + str(variant.end)
        if k in rec:
            variant.INFO["RUNREC"]=rec[k]

        w.write_record(variant)

    w.close()
    vcf.close()

