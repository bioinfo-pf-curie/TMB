#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of pyTMB software.
#
#  Copyright (c) 2020 - 2025 - Institut Curie
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
``pyTMB`` command-line entry point.

This module only handles argument parsing and wires CLI arguments to the
:func:`pytmb.tmb.calculate_tmb` library function.  All business logic lives
in :mod:`pytmb.tmb`.

Typical usage
-------------
After installation::

    pyTMB -i sample.vcf.gz \\
        --dbConfig config/snpeff.yml \\
        --varConfig config/mutect2.yml \\
        --effGenomeSize 33280000 \\
        --vaf 0.05 --maf 0.001 \\
        --minDepth 20 --minAltDepth 2 \\
        --filterLowQual --filterNonCoding --filterSyn --filterPolym
"""

import argparse
import sys

from pytmb import __version__
from pytmb.config import loadConfig
from pytmb.genome_size import getEffGenomeSizeFromBed
from pytmb.tmb import calculate_tmb, print_tmb_report


def argsParse():
    """
    Parse command-line arguments for the ``pyTMB`` command.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Calculate a Tumour Mutational Burden (TMB) score from a VCF file."
        ),
    )

    parser.add_argument(
        "-i", "--vcf", required=True,
        help="Input file (.vcf, .vcf.gz, .bcf, .bcf.gz)",
    )

    # Configs
    parser.add_argument(
        "--dbConfig", required=True,
        help="Databases config file (YAML)",
        type=str,
    )
    parser.add_argument(
        "--varConfig", required=True,
        help="Variant calling config file (YAML)",
        type=str,
    )
    parser.add_argument(
        "--sample",
        help="Specify the sample ID to focus on",
        type=str,
        default=None,
    )

    # Effective genome size
    parser.add_argument(
        "--effGenomeSize",
        help="Effective genome size (bp)",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--bed",
        help="Capture design BED file used to compute effGenomeSize when "
             "--effGenomeSize is not provided",
        default=None,
    )

    # Thresholds
    parser.add_argument(
        "--vaf",
        help="Filter variants with Allelic Ratio < vaf",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--maf",
        help="Filter variants with MAF > maf",
        type=float,
        default=1,
    )
    parser.add_argument(
        "--minDepth",
        help="Filter variants with depth < minDepth",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--minAltDepth",
        help="Filter variants with alternative allele depth < minAltDepth",
        type=int,
        default=1,
    )

    # Which variants to use
    parser.add_argument(
        "--filterLowQual",
        help="Filter low quality (i.e. not PASS) variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterIndels",
        help="Filter insertions/deletions",
        action="store_true",
    )
    parser.add_argument(
        "--filterCoding",
        help="Filter coding variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterSplice",
        help="Filter splice variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterNonCoding",
        help="Filter non-coding variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterSyn",
        help="Filter synonymous variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterNonSyn",
        help="Filter non-synonymous variants",
        action="store_true",
    )
    parser.add_argument(
        "--filterCancerHotspot",
        help="Filter variants annotated as cancer hotspots",
        action="store_true",
    )
    parser.add_argument(
        "--filterPolym",
        help="Filter polymorphism variants in genome databases (see --maf)",
        action="store_true",
    )
    parser.add_argument(
        "--filterRecurrence",
        help="Filter on run-level recurrence values",
        action="store_true",
    )

    # Databases
    parser.add_argument(
        "--polymDb",
        help="Databases used for polymorphism detection (comma-separated)",
        default="gnomad",
    )
    parser.add_argument(
        "--cancerDb",
        help="Databases used for cancer hotspot annotation (comma-separated)",
        default="cosmic",
    )

    # Misc
    parser.add_argument(
        "--verbose",
        help="Activate verbose mode",
        action="store_true",
    )
    parser.add_argument(
        "--debug",
        help="Export original VCF with TMB_FILTER tag",
        action="store_true",
    )
    parser.add_argument(
        "--export",
        help="Export a VCF with the passing variants to the specified path",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--version",
        help="Version number",
        action="version",
        version="%(prog)s (" + __version__ + ")",
    )

    return parser.parse_args()


def main():
    """
    Entry point for the ``pyTMB`` CLI command.
    """
    args = argsParse()

    # ------------------------------------------------------------------ #
    # Load configs                                                         #
    # ------------------------------------------------------------------ #
    db_flags = loadConfig(args.dbConfig)
    caller_flags = loadConfig(args.varConfig)

    # ------------------------------------------------------------------ #
    # Effective genome size                                                #
    # ------------------------------------------------------------------ #
    if args.effGenomeSize is None:
        if args.bed is not None:
            eff_gs = getEffGenomeSizeFromBed(args.bed, verbose=args.verbose)
        else:
            sys.stderr.write(
                "Error: Effective Genome Size not specified. "
                "See --effGenomeSize or --bed\n"
            )
            sys.exit(-1)
    else:
        eff_gs = args.effGenomeSize

    if eff_gs == 0:
        sys.stderr.write("Error: Effective genome size is 0. Check BED file.\n")
        sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Parameter validation                                                 #
    # ------------------------------------------------------------------ #
    if args.vaf > int(caller_flags["maxVaf"]):
        sys.stderr.write(
            "Error: vaf > maxVaf threshold. Check the configuration file.\n"
        )
        sys.exit(-1)

    if args.polymDb and not args.filterPolym:
        print(
            "Warning: --filterPolym argument not provided while --polymDb is specified!"
        )
    if args.polymDb and args.filterPolym and args.maf >= 1:
        print("Error: --maf must be < 1 !")
        sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Run TMB calculation                                                  #
    # ------------------------------------------------------------------ #
    results = calculate_tmb(
        vcf_path=args.vcf,
        db_flags=db_flags,
        caller_flags=caller_flags,
        eff_genome_size=eff_gs,
        sample=args.sample,
        vaf=args.vaf,
        maf=args.maf,
        min_depth=args.minDepth,
        min_alt_depth=args.minAltDepth,
        filter_low_qual=args.filterLowQual,
        filter_indels=args.filterIndels,
        filter_coding=args.filterCoding,
        filter_splice=args.filterSplice,
        filter_non_coding=args.filterNonCoding,
        filter_syn=args.filterSyn,
        filter_non_syn=args.filterNonSyn,
        filter_cancer_hotspot=args.filterCancerHotspot,
        filter_polym=args.filterPolym,
        filter_recurrence=args.filterRecurrence,
        polym_db=args.polymDb,
        cancer_db=args.cancerDb,
        export_path=args.export,
        debug=args.debug,
        verbose=args.verbose,
    )

    # ------------------------------------------------------------------ #
    # Print report                                                         #
    # ------------------------------------------------------------------ #
    print_tmb_report(results, args, __version__)


if __name__ == "__main__":
    main()
