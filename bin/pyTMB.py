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
Subset the annotation information to a few key values
"""
def subsetINFO(annot, keys):

    subsetInfo=[]
    for i in range(0, len(annotInfo)):
        z = dict((k, annotInfo[i][k]) for k in keys)
        subsetInfo.append(z)

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
Parse inputs
"""
def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")
    parser.add_argument("-c", "--config", help="Databases config file", default="./config/databases.yml")
    parser.add_argument("-a", "--annot", help="Annotation format used in the vcf file", default="snpeff")
    parser.add_argument("--bed", help="Capture design (BED file)", default=None)
    
    ## Thresholds
    parser.add_argument("--minVAF", help="Filter variants with Allelic Ratio < minAR", type=float, default=0.05)
    parser.add_argument("--minMAF", help="Filter variants with MAF < minMAF", type=float, default=0.001)
    parser.add_argument("--minDepth", help="Filter variants with depth < minDepth", type=int, default=5)
    
    ## Which variants to use
    parser.add_argument("--filterCoding", help="Filter Coding variants", action="store_true")
    parser.add_argument("--filterNonCoding", help="Filter Non-coding variants", action="store_true")
    parser.add_argument("--filterSyn", help="Filter Synonymous variants", action="store_true")
    parser.add_argument("--filterNonSyn", help="Filter Non-Synonymous variants", action="store_true")
    parser.add_argument("--filterHotspot", help="Filter variants in hotspots", action="store_true")
    parser.add_argument("--filterGermline", help="Filter potential variants flagged as germline in databases", action="store_true")
    
    ## Database
    parser.add_argument("--germlineDb", help="Databases used for germline annotation", default="gnomAD")
    parser.add_argument("--hotspotDb", help="Databases used for hotspot annotation", default="cosmic")
    
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
    
    for variant in VCF(args.vcf):
        varCounter+=1
        if (varCounter % 100 == 0):
            print ("## ",varCounter)
            if varCounter == 10000:
                sys.exit()

        try:
            ## Get annotInfo as a list of dict
            if args.annot == "snpeff":
                annotInfo = snpEff2dl(variant.INFO.get(dbflags['tag']))
            elif args.annot == "annovar" :
                annotInfo = annovar2dl(variant.INFO)

            ## No INFO field
            if annotInfo is None:
                continue

            if args.verbose:
                print (variant.CHROM, variant.start, variant.end, annotInfo)

            ## Variant Allele Frequency
            #if variant.INFO.get('AF') < args.minVAF:
            #    continue
                
            ## Sequencing Depth
            #if variant.INFO.get('DP') < args.minDepth:
            #    continue

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
            #if args.filterHotspot and isHotspot(v, infos=annotInfo, flags=dbflags['isHotspot'][args.hotspotDb]):
            #    continue

            ## Germline
            #if args.filterGermline and isGermline(v, infos=annotInfo, flags=dbflags['isGermline'][args.germlineDb]):
            #    continue

        except:
            ## TODO - try to catch different cases or at least to print info to check what's going on and where
            warnings.warn("Warning ...")
            raise

        ### Usage examples of cyvcf
        #variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']

        #variant.CHROM, variant.start, variant.end, variant.ID, \
        #    variant.FILTER, variant.QUAL
        
        # numpy arrays of specific things we pull from the sample fields.
        # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
        #variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
        #variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array
    
    
        ## INFO Field.
        ## extract from the info field by it's name:
        #variant.INFO.get('DP') # int
        #variant.INFO.get('FS') # float
        #variant.INFO.get('AC') # float
        
        # convert back to a string.
        #str(variant)
        
        ## per-sample info...
        # Get a numpy array of the depth per sample:
        #dp = variant.format('DP')
        # or of any other format field:
        #sb = variant.format('SB')
        #assert sb.shape == (n_samples, 4) # 4-values per

        ## Still alivie = use this variant for TMB calculation
        varTMB += 1

