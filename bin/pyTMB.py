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


def isCoding(v):
    return True

def isNonCoding(v):
    return True

def isSynonymous(v):
    return True

def isNonSynonymous(v):
    return True

def isHotspot(v):
    return True

def isSomatic(v):
    return True

def isGermline(v):
    return True


def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file (VCF format)")
    parser.add_argument("-b", "--bed", help="Capture design (BED file)", default='')
    
    ## Thresholds
    parser.add_argument("--minAR", help="Select variants with Allelic Ratio > minAR", default='')
    parser.add_argument("--minMAF", help="Select variants with MAF > minMAF", default='')
    parser.add_argument("--minDepth", help="Select variants with depth > minDepth", default=5)
    
    ## Which variants to use ?
    parser.add_argument("--useCoding", help="Include Coding variants", default='')
    parser.add_argument("--useNonCoding", help="Include Non-coding variants", default='')
    parser.add_argument("--useSyn", help="Include Synonymous variants", default='')
    parser.add_argument("--useNonSyn", help="Include Non-Synonymous variants", default='')
    parser.add_argument("--useHotspot", help="Include variants in hotspots", default='')
    parser.add_argument("--filterGermline", help="Filter potential variants flagged as germline in databases", default='')
    
    ## Database
    parser.add_argument("--germlineDb", help="Databases used for germline annotation", default='')
    parser.add_argument("--hotspotDb", help="Databases used for hotspot annotation", default='')
    
    ## Others
    parser.add_argument("-v", "--version", help="Pipeline version", default='')
    
    args = parser.parse_args()
    return (args)



if __name__ == "__main__":

    args = argsParse()

    for variant in VCF(args.input): # or VCF('some.bcf')

        if variant.INFO.get('AF') > args.minMAF:
            continue

        if variant.INFO.get('DP') > args.minDepth:
            continue


        ### Usage examples of cyvcf
        variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']

        variant.CHROM, variant.start, variant.end, variant.ID, \
            variant.FILTER, variant.QUAL
        
        # numpy arrays of specific things we pull from the sample fields.
        # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
        variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
        variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array
    
    
        ## INFO Field.
        ## extract from the info field by it's name:
        variant.INFO.get('DP') # int
        variant.INFO.get('FS') # float
        variant.INFO.get('AC') # float
        
        # convert back to a string.
        str(variant)
        
        ## per-sample info...
        # Get a numpy array of the depth per sample:
        dp = variant.format('DP')
        # or of any other format field:
        sb = variant.format('SB')
        assert sb.shape == (n_samples, 4) # 4-values per
        
