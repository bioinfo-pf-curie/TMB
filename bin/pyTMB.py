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

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file (VCF format)")
parser.add_argument("-b", "--bed", help="Capture design (BED file)", default='')

## Thresholds
parser.add_argument("-r", "--minAR", help="", default='')
parser.add_argument("-m", "--minMAF", help="", default='') 
parser.add_argument("-c", "--minCov", help="", default='') 

## Which variants to use ?
parser.add_argument("-uc", "--useCoding", help="", default='')
parser.add_argument("-unc", "--useNonCoding", help="", default='') 
parser.add_argument("-us", "--useSyn", help="", default='') 
parser.add_argument("-uns", "--useNonSyn", help="", default='') 
parser.add_argument("-us", "--useSomatic", help="", default='')
parser.add_argument("-ug", "--useGermline", help="", default='') 
parser.add_argument("-uh", "--useHotspot", help="", default='') 

## Database
parser.add_argument("-ds", "--somaticDb", help="", default='') 
parser.add_argument("-dh", "--hotspotDb", help="", default='') 

## Others
parser.add_argument("-v", "--version", help="Pipeline version", default='')


args = parser.parse_args()


for variant in VCF(args.input): # or VCF('some.bcf')
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

