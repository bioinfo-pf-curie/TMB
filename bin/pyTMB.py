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



"""
Step of the script:  
- Read and parse annotated vcf file (PySam) 
- TAG each variant directly in the .vcf file  

- Variants with Allelic Ratio < 5%: LowFreq/PASS ## --minAR
- Variants with MAF < XX: LowFreq/PASS ## --minMAF
- Variants with Coverage < XX: LowDepth/PASS ## --minCov

- Variants annotated as coding/non coding: Relevant /NonRelevant ## --useNonCoding / --useCoding
- Synonymous variant: NonSyn/Syn ## --useSyn / --useNonSyn 
- Variants found in Gnomad: Germline/Somatic ## --useSomatic / --useGermline / --somaticDb
- Variants found in COSMIC: Onco/PASS ## --useHotspot / --hotspotDb


- Select variant of interest only Relevant, NonSyn, Somatic and PASS tagged variants and filter other variants (NonRelevant, Syn, Germline, LowFreq, LowDepth...) 
- Extract exome size (bp) 
- Count total number of variants 
- Compute the TMB score  
"""


from cyvcf2 import VCF
import argparse
import sys
import warnings
import re
import yaml


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
            sys.stderr.write("Error : wrong input format in line", nline,". Not a BED file !?")
            sys.exit(-1)

        intl = abs(int(end) - int(start))
        effgs += intl
    
    bed_handle.close()
    return effgs


"""
Check if a variant has the provided annotation flags
"""
def isAnnotatedAs(v, infos, flags):

    ## Subset annotation information as list
    subINFO = subsetINFO(infos, keys=flags.keys())
    
    ## Compare list variant information and expected list of tags
    for sub in subINFO:
        ## For all keys
        for k in flags.keys():
            ## For each value in keys
            for val in flags[k]:
                ## Use search to deal with multiple annotations
                if re.search(val, sub[k]):
                    return(True)

    return(False)


"""
Check if a variant is in a genomeDb with a MAF > val
"""
def isGenomeDb(v, infos, flags, val):

    subINFO = subsetINFO(infos, keys=flags)
    for key in subINFO:
        if subINFO[key] is not None and subINFO[key] != ".":
            if float(subINFO[key]) > float(val):
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
        subsetInfo=[]
        for i in range(0, len(annot)):
            z = dict((k, annot[i][k]) for k in keys if k in annot)
            subsetInfo.append(z)
    else:
        subsetInfo = dict((k, annot[k]) for k in keys if k in annot)

    return(subsetInfo)

"""
Format the INFO field from snpEff and return a list of dict
"""
def snpEff2dl(INFO):

    if INFO is not None:
        annotTag = INFO.split(',')                                                                                                                                       
        annotInfo=[]                                                                                                                                                                                 
        for i in range(0, len(annotTag)):                                                                                                                                                            
            annot=annotTag[i].split('|')                                                                                                                                                             
            dictannot = {i : annot[i] for i in range(0, len(annot))}                                                                                                                                 
            annotInfo.append(dictannot)                                                                                                                                                              
        return(annotInfo)

"""
Format the INFO field from ANNOVAR and return a list of dict
"""
def annovar2dl(INFO):

    if INFO is not None:
        return [dict(INFO)]

"""
Get a tag value from either the format field or the info field
"""
def getTag(v, tag):

    ## First check in FORMAT field
    if tag in variant.FORMAT:
        val=variant.format(tag, str)
    ## Otherwise, check in INFO field
    if tag not in variant.FORMAT or val is None:
        val=variant.INFO.get(tag)

    return float(val)

"""
Parse inputs
"""
def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")
    parser.add_argument("-c", "--config", help="Databases config file", default="./config/databases.yml")
    parser.add_argument("-a", "--annot", help="Annotation format used in the vcf file", default="snpeff")

    ## Efective genome size
    parser.add_argument("--effGenomeSize", help="Effective genome size", type=int, default=None)
    parser.add_argument("--bed", help="Capture design to use if effGenomeSize is not defined (BED file)", default=None)
    
    ## Thresholds
    parser.add_argument("--minVAF", help="Filter variants with Allelic Ratio < minAR", type=float, default=0.05)
    parser.add_argument("--minMAF", help="Filter variants with MAF < minMAF", type=float, default=0.001)
    parser.add_argument("--minDepth", help="Filter variants with depth < minDepth", type=int, default=5)
    
    ## Which variants to use
    parser.add_argument("--filterLowQual", help="Filter low quality (i.e not PASS) variant", action="store_true")
    parser.add_argument("--filterCoding", help="Filter Coding variants", action="store_true")
    parser.add_argument("--filterNonCoding", help="Filter Non-coding variants", action="store_true")
    parser.add_argument("--filterSyn", help="Filter Synonymous variants", action="store_true")
    parser.add_argument("--filterNonSyn", help="Filter Non-Synonymous variants", action="store_true")
    parser.add_argument("--filterCancerHotspot", help="Filter variants annotated as cancer hotspots", action="store_true")
    parser.add_argument("--filterGenomeDb", help="Filter variants detected in large genome databases. See --minMAF", action="store_true")
    
    ## Databases
    parser.add_argument("--genomeDb", help="Databases used for large genome annotation (comma separated)", default="gnomad")
    parser.add_argument("--cancerDb", help="Databases used for cancer hotspot annotation", default="cosmic")
    
    ## Others
    parser.add_argument("-v", "--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--version", help="Pipeline version", action="store_true")
    
    args = parser.parse_args()
    return (args)


if __name__ == "__main__":

    args = argsParse()
    dbconfig = loadConfig(args.config)
    dbflags = dbconfig[args.annot]

    varCounter = 0
    varTMB = 0

    if args.effGenomeSize is None:
        if args.bed is not None:
            effGS = getEffGenomeSizeFromBed(args.bed)
        else:
            sys.stderr.write("Error: Effective Genome Size not specified. See --effGenomeSize or --bed")
            sys.exit(-1)
    else:
        effGS = args.effGenomeSize
    

    for variant in VCF(args.vcf):
        varCounter+=1
        if (varCounter % 1000 == 0):
            print ("## ",varCounter)
            if varCounter == 100000:
                sys.exit()

        try:
            ## All vcf INFO
            dbInfo=dict(variant.INFO)

            ## Get annotation INFO as a list of dict
            if args.annot == "snpeff":
                annotInfo = snpEff2dl(variant.INFO.get(dbflags['annotTag']))
            elif args.annot == "annovar" :
                annotInfo = annovar2dl(variant.INFO)

            ## No INFO field
            if dbInfo is None or annotInfo is None:
                continue
            
            if args.verbose:
                print (variant.CHROM, variant.start, variant.end, variant.QUAL, variant.FILTER, dbInfo, variant.FORMAT)

            ## Variant Allele Frequency
            if getTag(variant, dbflags['freqTag']) < args.minVAF:
                print("FILTER FREQ")
                continue
                
            ## Sequencing Depth
            if getTag(variant, dbflags['depthTag']) < args.minDepth:
                print("FILTER DEPTH")
                continue

            ## Variant has a QUAL value or not PASS in the FILTER column
            if args.filterLowQual and (variant.QUAL is not None or variant.FILTER is not None):
                continue
 
            ## Coding variants
            if args.filterCoding and isAnnotatedAs(variant, infos=annotInfo, flags=dbflags['isCoding']):
                print("FILTER IS_CODING")
                continue

            ## Non-coding variants
            if args.filterNonCoding and isAnnotatedAs(variant, infos=annotInfo, flags=dbflags['isNonCoding']):
                print("FILTER IS_NON_CODING")
                continue 

            ## Synonymous
            if args.filterSyn and isAnnotatedAs(variant, infos=annotInfo, flags=dbflags['isSynonymous']):
                print("FILTER IS_SYNONYMOUS")
                continue
            
            ## Non synonymous
            if args.filterNonSyn and isAnnotatedAs(variant, infos=annotInfo, flags=dbflags['isNonSynonymous']):
                print("FILTER IS_NON_SYNONYMOUS")
                continue

            ## Hotpost
            if args.filterCancerHotspot:
                ## Flatten list of fields
                fdb=[]
                for db in args.cancerDb.split(','):
                    for x in dbflags['cancerDb'][db]:
                        fdb.append(x)
 
                if isCancerHotspot(variant, infos=dbInfo, flags=fdb):
                    print("FILTER HOSTPOT")
                    continue

            ## Large Genome Cohort
            if args.filterGenomeDb :
                ## Flatten list of fields
                fdb=[]
                for db in args.genomeDb.split(','):
                    for x in dbflags['genomeDb'][db]:
                        fdb.append(x)

                if isGenomeDb(variant, infos=dbInfo, flags=fdb, val=args.minMAF):
                    print("FILTER MAF")
                    continue

                ## Still alive
                varTMB += 1
        except:
            ## TODO - try to catch different cases or at least to print info to check what's going on and where
            warnings.warn("Warning ...")
            raise

        ## Calculate TMB
        TMB = float(varTMB)/float(effGS)*1e6
        print("TMB=",TMB)

