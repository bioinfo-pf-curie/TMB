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


"""
Check if the variant affects the coding sequence
"""
def isCoding(v, db="snpEff"):
    ret=False

    ## snpEff annotation
    snpEffFlag=["intragenic_variant", "5_prime_UTR", "3_prime_UTR"]
    ## ANNOVAR annotation
    annovarFlag=[]

    if db == "snpEff":
        flags = snpEffFlag
    elif db == "annovar":
        flags = annovarFlag

    annotInfo=v.INFO.get('ANN').split(',')
    for i in range(0, len(annotInfo)):
        annot=annotInfo[i].split('|')
        for pattern in flags:
            p=re.compile(pattern)
            if p.match(annot[1]):
                ret=True
                break
        
    if ret:
        print("CODING")

    return ret


"""
Check if the variant is within a non-coding region
"""
def isNonCoding(v):
    ret=False
    
    ## snpEff annotation
    flags=["intergenic_region", "intron_variant", "upstream_gene_variant"]

    annotInfo=v.INFO.get('ANN').split(',')
    for i in range(0, len(annotInfo)):
        annot=annotInfo[i].split('|')
        for pattern in flags:
            p=re.compile(pattern)
            if p.match(annot[1]):
                ret=True
                break

    if ret:
        print("NONCODING")

    return ret


"""
Is synonymous variant
"""
def isSynonymous(v):
    ret=False

    ## snpEff annotation
    flags=["synonymous_variant"]

    annotInfo=v.INFO.get('ANN').split(',')
    for i in range(0, len(annotInfo)):
        annot=annotInfo[i].split('|')
        for pattern in flags:
            p=re.compile(pattern)
            if p.match(annot[1]):
                ret=True
                break
    if ret:
        print("SYN")

    return ret

"""
Is mis-sense variants
"""
def isNonSynonymous(v):
    ret=False

    ## snpEff annotation
    flags=["missense_variant", "initiator_codon_variant", "stop_retained_variant", "stop_gained"]

    annotInfo=v.INFO.get('ANN').split(',')
    for i in range(0, len(annotInfo)):
        annot=annotInfo[i].split('|')
        for pattern in flags:
            p=re.compile(pattern)
            if p.match(annot[1]):
                ret=True
                break
    if ret:
        print("NONSYN")

    return ret


def isHotspot(v, databases):
    print (v.INFO.get('cosmic_coding_ID'))
    print (v.INFO.get('cosmic_noncoding_ID'))

    return True

def isGermline(v, databases):
    return True


def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")

    parser.add_argument("--bed", help="Capture design (BED file)", default=None)
    
    ## Thresholds
    parser.add_argument("--minVAF", help="Select variants with Allelic Ratio > minAR", type=float, default=0.05)
    parser.add_argument("--minMAF", help="Select variants with MAF > minMAF", type=float, default=0.001)
    parser.add_argument("--minDepth", help="Select variants with depth > minDepth", type=int, default=5)
    
    ## Which variants to use
    parser.add_argument("--useCoding", help="Include Coding variants", action="store_true")
    parser.add_argument("--useNonCoding", help="Include Non-coding variants", action="store_true")
    parser.add_argument("--useSyn", help="Include Synonymous variants", action="store_true")
    parser.add_argument("--useNonSyn", help="Include Non-Synonymous variants", action="store_true")

    ## Which variants to filter
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
    varCounter=0
    varTMB=0
    
    for variant in VCF(args.vcf):
        varCounter+=1
        if (varCounter % 100 == 0):
            #print ("## ",varCounter)
            if varCounter == 100000:
                sys.exit()

        if args.verbose:
            print (variant.CHROM, variant.start, variant.end, variant.INFO.get("ANN"))

        try:
            ## Variant Allele Frequency
            if variant.INFO.get('AF') < args.minVAF:
                continue
                
            ## Sequencing Depth
            if variant.INFO.get('DP') < args.minDepth:
                continue
                    
            ## MAF
                    
            ## Coding variants
            if args.useCoding and not isCoding(variant):
                continue

            ## Non-coding variants
            if args.useNonCoding and not isNonCoding(variant):
                continue 

            ## Synonymous
            if args.useSyn and not isSynonymous(variant):
                continue
            
            ## Non synonymous
            if args.useNonSyn and not isNonSynonymous(v):
                continue

            ## Hotpost
            if args.filterHotspot and isHotspot(v, args.hotspotDb):
                continue

            ## Germline
            if args.filterGermline and isGermline(v, args.germlineDb):
                continue
                
        except:
            ## TODO - try to catch different cases or at least to print info to check what's going on and where
            warnings.warn("Warning ...")
        

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

