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

__version__ = '1.3.0dev'

"""
This script is designed to calculate a TMB score from a VCF file.
If no filters are specified, all variants will be used.

Let's defined a TMB as a score using PASS, non-synonymous, coding, non polymorphism variants ...
In this case, a typical usage would be :

python pyTMB.py -i ${VCF} --effGenomeSize 33280000 \
--vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2\
--filterLowQual \
--filterNonCoding \
--filterSyn \
--filterPolym \
--polymDb 1k,gnomad \
--dbConfig ${DB_CONFIG} \
--varConfig ${VAR_CONFIG} > TMB_results.log
"""

import cyvcf2
import argparse
import sys
import warnings
import re
import yaml
import numpy as np
import os.path
from datetime import date

"""
Load yaml file
"""


def loadConfig(infile):

    with open(infile, 'r') as stream:
        try:
            return(yaml.safe_load(stream))
        except:
            raise


"""
Extract VCF column name with multi-Allelic information
"""


def getMultiAlleleHeader(vcf):

    FORMAT = []
    INFO = []
    for h in vcf.header_iter():
        i = h.info(extra=True)
        if 'Number' in i.keys() and i['Number'] == 'A':
            if i['HeaderType'] == 'FORMAT':
                FORMAT.append(i['ID'])
            elif i['HeaderType'] == 'INFO':
                INFO.append(i['ID'])
    return(dict(FORMAT=FORMAT, INFO=INFO))


"""
Calculate Effective Genome Size from a BED file
"""


def getEffGenomeSizeFromBed(infile, verbose=False):

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


"""
Check if a variant has the provided annotation flags
"""


def isAnnotatedAs(v, infos, flags, sep):

    # Subset annotation information as list
    subINFO = subsetINFO(infos, keys=flags.keys())
    # Compare list variant information and expected list of tags
    for sub in subINFO:
        # For all keys
        for k in flags.keys():
            # For each value in keys
            for val in flags[k]:
                # Use search to deal with multiple annotations
                for subval in sub[k].split(sep):
                    if val == subval:
                        return(True)
    return(False)


"""
Check if a variant is in a genomeDb with a MAF > val
"""


def isPolym(v, infos, flags, val):

    subINFO = subsetINFO(infos, keys=flags)
    for key in subINFO:
        if type(subINFO[key]) is tuple:
            for i in subINFO[key]:
                if i is not None and i != ".":
                    if float(i) >= float(val):
                        return True
        elif subINFO[key] is not None and subINFO[key] != ".":
            if float(subINFO[key]) >= float(val):
                return True
    return(False)


"""
Check if a variant is annotated as a cancer hotspot
"""


def isCancerHotspot(v, infos, flags):

    subINFO = subsetINFO(infos, keys=flags)
    for key in subINFO:
        if subINFO[key] is not None and subINFO[key] != ".":
            return True
    return(False)


"""
Subset the annotation information to a few key values
"""


def subsetINFO(annot, keys):

    if isinstance(annot, list):
        subsetInfo = []
        for i in range(0, len(annot)):
            z = dict((k, annot[i][k]) for k in keys if k in annot[i])
            if len(z) > 0:
                subsetInfo.append(z)
    else:
        subsetInfo = dict((k, annot[k]) for k in keys if k in annot)
    return(subsetInfo)


"""
Format the INFO field from snpEff and return a list of dict
ie. snpEff
"""


def infoTag2dl(INFO):

    if INFO is not None:
        annotTag = INFO.split(',')
        annotInfo = []
        for i in range(0, len(annotTag)):
            annot = annotTag[i].split('|')
            dictannot = {i: annot[i] for i in range(0, len(annot))}
            annotInfo.append(dictannot)
        return(annotInfo)


"""
Format the INFO field from ANNOVAR and return a list of dict
ie. annovar
"""


def info2dl(INFO):

    if INFO is not None:
        return [dict(INFO)]


"""
Get a tag value from either the format field or the info field
Return a 2D numpy array
"""


def getTag(v, tag):

    # First check in FORMAT field
    if tag in variant.FORMAT:
        val = variant.format(tag)
        # if np.shape(val) == (1, 1) and val < 0:
        #        val = None

    # Otherwise, check in INFO field
    if tag not in variant.FORMAT or val is None:
        val = variant.INFO.get(tag)

    # if np.shape(val) == (1, 1) and val is not None:
    #    val = float(val)

    if type(val) != np.ndarray:
        val = np.array([val], float)
    else:
        val = val.astype('float')

    return(val)


"""
Parse inputs
"""


def argsParse():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")

    # Configs
    parser.add_argument("--dbConfig", help="Databases config file", type=str, default="./config/databases.yml")
    parser.add_argument("--varConfig", help="Variant calling config file", type=str, default="./config/calling.yml")
    parser.add_argument("--sample", help="Specify the sample ID to focus on", type=str, default=None)

    # Efective genome size
    parser.add_argument("--effGenomeSize", help="Effective genome size", type=int, default=None)
    parser.add_argument("--bed", help="Capture design to use if effGenomeSize is not defined (BED file)", default=None)

    # Thresholds
    parser.add_argument("--vaf", help="Filter variants with Allelic Ratio <= vaf", type=float, default=0.05)
    parser.add_argument("--maf", help="Filter variants with MAF > maf", type=float, default=0.001)
    parser.add_argument("--minDepth", help="Filter variants with depth < minDepth", type=int, default=5)
    parser.add_argument("--minAltDepth", help="Filter variants with alternative allele depth <= minAltDepth", type=int, default=2)

    # Which variants to use
    parser.add_argument("--filterLowQual", help="Filter low quality (i.e not PASS) variant", action="store_true")
    parser.add_argument("--filterIndels", help="Filter insertions/deletions", action="store_true")
    parser.add_argument("--filterCoding", help="Filter Coding variants", action="store_true")
    parser.add_argument("--filterSplice", help="Filter Splice variants", action="store_true")
    parser.add_argument("--filterNonCoding", help="Filter Non-coding variants", action="store_true")
    parser.add_argument("--filterSyn", help="Filter Synonymous variants", action="store_true")
    parser.add_argument("--filterNonSyn", help="Filter Non-Synonymous variants", action="store_true")
    parser.add_argument("--filterCancerHotspot", help="Filter variants annotated as cancer hotspots", action="store_true")
    parser.add_argument("--filterPolym", help="Filter polymorphism variants in genome databases. See --maf", action="store_true")
    parser.add_argument("--filterRecurrence", help="Filter on recurrence values", action="store_true")

    # Databases
    parser.add_argument("--polymDb", help="Databases used for polymorphisms detection (comma separated)", default="gnomad")
    parser.add_argument("--cancerDb", help="Databases used for cancer hotspot annotation (comma separated)", default="cosmic")

    # Others
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--debug", help="Export original VCF with TMB_FILTER tag", action="store_true")
    parser.add_argument("--export", help="Export a VCF with the considered variants", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s ("+__version__+")")

    args = parser.parse_args()
    return (args)


if __name__ == "__main__":

    args = argsParse()

    # Load Data
    if args.sample is not None:
        vcf = cyvcf2.VCF(args.vcf)
        count = 0
        for sample in vcf.samples:
            count = count + 1
            if str(sample) == str(args.sample):
                vcf = cyvcf2.VCF(args.vcf, samples = args.sample)
            elif count == len(vcf.samples):
                print("Error: Name of the sample incorrect")
                sys.exit(-1)
    else:
        vcf = cyvcf2.VCF(args.vcf)

    # Sample name
    if len(vcf.samples) > 1:
        sys.stderr.write("Error: " + str(len(vcf.samples)) +
                         " sample detected. This version is designed for a single sample ! \n")
        sys.exit(-1)

    # Ouptuts
    if args.export:
        wx = cyvcf2.Writer(re.sub(r'\.vcf$|\.vcf.gz$|\.bcf',
                                  '_export.vcf', os.path.basename(args.vcf)), vcf)
    if args.debug:
        vcf.add_info_to_header({'ID': 'TMB_FILTERS', 'Description': 'Detected filters for TMB calculation',
                                'Type': 'Character', 'Number': '1'})
        wd = cyvcf2.Writer(re.sub(r'\.vcf$|\.vcf.gz$|\.bcf',
                                  '_debug.vcf', os.path.basename(args.vcf)), vcf)

    # Load config
    dbFlags = loadConfig(args.dbConfig)
    callerFlags = loadConfig(args.varConfig)

    # Genome size
    if args.effGenomeSize is None:
        if args.bed is not None:
            effGS = getEffGenomeSizeFromBed(args.bed)
        else:
            sys.stderr.write(
                "Error: Effective Genome Size not specified. See --effGenomeSize or --bed")
            sys.exit(-1)
    else:
        effGS = args.effGenomeSize


    varCounter = 0
    varNI = 0
    varTMB = 0
    for variant in vcf:
        varCounter += 1
        if (varCounter % 1000 == 0 and args.verbose):
            print ("## ", varCounter)
            if args.debug and varCounter == 1000:
                sys.exit()

        try:
            # All vcf INFO
            dbInfo = dict(variant.INFO)
            debugInfo = ""

            # Get annotation INFO as a list of dict
            if dbFlags['tag'] != '':
                annotInfo = infoTag2dl(variant.INFO.get(dbFlags['tag']))

            else:
                annotInfo = info2dl(variant.INFO)

            # No INFO field
            if dbInfo is None or annotInfo is None:
                varNI += 1
                continue

            # Indels
            if args.filterIndels and variant.is_indel:
                debugInfo = ",".join([debugInfo, "INDEL"])
                if not args.debug:
                    continue

            # Variant has a QUAL value or not PASS in the FILTER column
            if args.filterLowQual and (variant.QUAL is not None or variant.FILTER is not None):
                debugInfo = ",".join([debugInfo, "QUAL"])
                if not args.debug:
                    continue

            #######################
            # FORMAT
            #######################

            # Variant Allele Frequency
            fval = getTag(variant, callerFlags['freq'])
            if fval is not None and len(fval[fval <= args.vaf]) == len(variant.ALT):
                debugInfo = ",".join([debugInfo, "VAF"])
                if not args.debug:
                    continue

            # Sequencing Depth
            dval = getTag(variant, callerFlags['depth'])
            if dval is not None and len(dval[dval <= args.minDepth]) == len(variant.ALT):
                debugInfo = ",".join([debugInfo, "DEPTH"])
                if not args.debug:
                    continue

            # Alternative allele Depth
            ad = getTag(variant, callerFlags['altDepth'])
            # case where AD = REF + ALTs
            if len(ad) == (len(variant.ALT) + 1):
                ad = ad[1:]

            if ad is not None and len(ad[ad <= args.minAltDepth]) == len(variant.ALT):
                debugInfo = ",".join([debugInfo, "ALTDEPTH"])
                if not args.debug:
                    continue

            ######################
            # Annotation INFO
            ######################

            # Coding variants
            if args.filterCoding and isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isCoding'], sep=dbFlags['sep']):
                if not isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isNonCoding'], sep=dbFlags['sep']):
                    debugInfo = ",".join([debugInfo, "CODING"])
                    if not args.debug:
                        continue

            # Splice variants
            if args.filterSplice and isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isSplicing'], sep=dbFlags['sep']):
                debugInfo = ",".join([debugInfo, "SPLICING"])
                if not args.debug:
                    continue

            # Non-coding variants
            if args.filterNonCoding and isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isNonCoding'], sep=dbFlags['sep']):
                 if not isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isCoding'], sep=dbFlags['sep']):
                    debugInfo = ",".join([debugInfo, "NONCODING"])
                    if not args.debug:
                        continue

            # Synonymous
            if args.filterSyn and isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isSynonymous'], sep=dbFlags['sep']):
                if not isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isNonSynonymous'], sep=dbFlags['sep']):
                    debugInfo = ",".join([debugInfo, "SYN"])
                    if not args.debug:
                        continue

            # Non synonymous
            if args.filterNonSyn and isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isNonSynonymous'], sep=dbFlags['sep']):
                if not isAnnotatedAs(variant, infos=annotInfo, flags=dbFlags['isSynonymous'], sep=dbFlags['sep']):
                    debugInfo = ",".join([debugInfo, "NON_SYN"])
                    if not args.debug:
                        continue

            # Hotspot
            if args.filterCancerHotspot:
                # Flatten list of fields
                fdb = []
                for db in args.cancerDb.split(','):
                    for x in dbFlags['cancerDb'][db]:
                        fdb.append(x)

                if isCancerHotspot(variant, infos=dbInfo, flags=fdb):
                    debugInfo = ",".join([debugInfo, "HOTSPOT"])
                    if not args.debug:
                        continue

            # Polymorphisms
            if args.filterPolym:
                # Flatten list of fields
                fdb = []
                for db in args.polymDb.split(','):
                    for x in dbFlags['polymDb'][db]:
                        fdb.append(x)

                if isPolym(variant, infos=dbInfo, flags=fdb, val=args.maf):
                    debugInfo = ",".join([debugInfo, "POLYM"])
                    if not args.debug:
                        continue

            # Recurrence
            if args.filterRecurrence:
                if isPolym(variant, infos=dbInfo, flags=dbFlags['recurrence']['run'], val=0):
                    debugInfo = ",".join([debugInfo, "RUNREC"])
                    if not args.debug:
                        continue

        except:
            warnflag = str(variant.CHROM) + ":" + str(variant.start) + "-" + str(variant.end)
            warnings.warn("Warning : variant ", warnflag, " raises an error. Skipped so far ...")
            raise

        # Still alive
        if debugInfo == "":
            varTMB += 1
            if args.export:
                wx.write_record(variant)

        if args.debug:
            variant.INFO["TMB_FILTERS"] = re.sub(r'^,', '', debugInfo)
            wd.write_record(variant)

    if args.export:
        wx.close()
    if args.debug:
        wd.close()
    vcf.close()

    # Calculate TMB
    TMB = round(float(varTMB)/(float(effGS)/1e6), 2)

    # No Annotation
    if varNI == varCounter:
        print("/!\ Warning: No Annotation detected. Is your file annotated ?")
        sys.exit(-1)

    # Output
    print("pyTMB version=", __version__)
    print("When=", date.today())
    print("")
    print("Config caller=", args.varConfig)
    print("Config databases=", args.dbConfig)
    print("")
    print("Filters:")
    print("-------")
    print("VAF=", args.vaf)
    print("MAF=", args.maf)
    print("minDepth=", args.minDepth)
    print("minAltDepth=", args.minAltDepth)
    print("filterLowQual=", args.filterLowQual)
    print("filterIndels=", args.filterIndels)
    print("filterCoding=", args.filterCoding)
    print("filterNonCoding=", args.filterNonCoding)
    print("filterSplice=", args.filterSplice)
    print("filterSyn=", args.filterSyn)
    print("filterNonSyn=", args.filterNonSyn)
    print("filterConcerHostpot=", args.filterCancerHotspot)
    print("filterPolym=", args.filterPolym)
    print("")
    print("Total number of variants=", varCounter)
    print("Non-informative variants=", varNI)
    print("Variants after filters=", varTMB)
    print("Effective Genome Size=", effGS)
    print("")
    print("TMB=", TMB)
