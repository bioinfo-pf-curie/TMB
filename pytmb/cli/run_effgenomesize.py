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
``pyEffGenomeSize`` command-line entry point.

This script calculates an effective genome size from:

- A capture BED file
- A genome annotation GTF file
- (optionally) a BAM file via mosdepth for coverage-based filtering

Typical usage::

    pyEffGenomeSize --bed design.bed --gtf genome.gtf \\
        --filterNonCoding --oprefix results/effgs
"""

import argparse
import os
import sys
import subprocess as sbp

from pytmb import __version__
from pytmb.genome_size import getEffGenomeSizeFromMosdepth


# ---------------------------------------------------------------------------
# GTF filtering helpers (pybedtools-based)
# ---------------------------------------------------------------------------

def filterGtf(interval, featuretype):
    """
    pybedtools filter callback: keep intervals whose ``transcript_type``
    attribute is in *featuretype*.

    Parameters
    ----------
    interval : pybedtools.Interval
    featuretype : list of str

    Returns
    -------
    str or False
    """
    if "transcript_type" in interval.attrs:
        if interval.attrs["transcript_type"] in featuretype:
            return interval.attrs["transcript_type"]
    return False


def filterFeatureGtf(interval, features):
    """
    pybedtools filter callback: keep intervals whose feature type (3rd column)
    is in *features*.

    Parameters
    ----------
    interval : pybedtools.Interval
    features : list of str

    Returns
    -------
    str or False
    """
    if interval[2] in features:
        return interval[2]
    return False


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def argsParse():
    """
    Parse command-line arguments for the ``pyEffGenomeSize`` command.

    Returns
    -------
    argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Calculate an effective genome size from a capture BED file, "
            "a GTF annotation, and optionally a BAM file via mosdepth."
        ),
    )

    # Inputs
    parser.add_argument(
        "--bed", help="BED file (.bed)", required=True, default=None
    )
    parser.add_argument(
        "--gtf", help="GTF file for genome annotation (.gtf)", required=True, default=None
    )
    parser.add_argument(
        "--bam", help="BAM file for mapping statistics (.bam)", default=None
    )

    # Filters
    parser.add_argument(
        "--mosdepth",
        help="Enable mosdepth processing (requires --bam)",
        action="store_true",
    )
    parser.add_argument(
        "--minCoverage",
        help="Minimum coverage of the region",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--minMapq",
        help="Minimum mapping quality",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--filterNonCoding",
        help="Filter regions associated with non-coding annotations",
        action="store_true",
    )
    parser.add_argument(
        "--filterCoding",
        help="Filter regions associated with coding annotations",
        action="store_true",
    )
    parser.add_argument(
        "--featureTypes",
        help=(
            "List of feature types (exon, gene, transcript, UTR, CDS) to keep "
            "(3rd GTF column).  Required with --filterCoding."
        ),
        nargs="+",
        default=[],
    )

    # Output / misc
    parser.add_argument(
        "--saveIntermediates",
        help="Keep mosdepth intermediate files",
        action="store_true",
    )
    parser.add_argument(
        "-t", "--thread",
        help="Number of threads for mosdepth",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--oprefix",
        help="Output file prefix",
        type=str,
        default="pyeffg",
    )
    parser.add_argument(
        "--verbose",
        help="Activate verbose mode",
        action="store_true",
    )
    parser.add_argument(
        "--version",
        help="Version number",
        action="version",
        version="%(prog)s (" + __version__ + ")",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """
    Entry point for the ``pyEffGenomeSize`` CLI command.
    """
    args = argsParse()

    # Lazy import: pybedtools is optional / heavy
    try:
        import pybedtools as pbt
        import pandas as pd
    except ImportError as exc:
        sys.stderr.write(
            f"Error: missing dependency for pyEffGenomeSize: {exc}\n"
            "Install pybedtools and pandas to use this command.\n"
        )
        sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Validations                                                          #
    # ------------------------------------------------------------------ #
    if (args.bam or args.minCoverage or args.minMapq) and not args.mosdepth:
        sys.stderr.write(
            "Error: --mosdepth is required if --bam, --minCoverage or --minMapq "
            "is used.\n"
        )
        sys.exit(-1)

    if not (args.filterNonCoding or args.filterCoding):
        sys.stderr.write(
            "Error: --filterCoding or --filterNonCoding is required!\n"
        )
        sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Load BED / GTF                                                       #
    # ------------------------------------------------------------------ #
    if args.verbose:
        print("[RUNNING INFO]: Loading data\n")

    my_bed = pbt.BedTool(args.bed)
    my_gtf = pbt.BedTool(args.gtf)

    if str(my_gtf[2]).find("transcript_type") == -1:
        sys.stderr.write(
            "Error: GTF doesn't have transcript_type info! Can't filter this file.\n"
        )
        sys.exit(-1)

    # Warn about annotations missing transcript_type
    count = sum(1 for interval in my_gtf if "transcript_type" not in interval.attrs)
    if count > 0:
        print(
            f"Warning: {count} annotations not used because no transcript_type info"
        )

    # ------------------------------------------------------------------ #
    # Define coding / non-coding transcript type sets                     #
    # ------------------------------------------------------------------ #
    coding = [
        "antisense", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
        "protein_coding", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
        "nonsense_mediated_decay", "non_stop_decay",
    ]
    non_coding = [
        "3prime_overlapping_ncrna", "IG_C_pseudogene", "IG_J_pseudogene",
        "IG_V_pseudogene", "lincRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA",
        "polymorphic_pseudogene", "processed_transcript", "processed_pseudogene",
        "pseudogene", "retained_intron", "rRNA", "sense_intronic",
        "sense_overlapping", "snoRNA", "snRNA",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
        "unitary_pseudogene", "unprocessed_pseudogene",
    ]

    if args.verbose:
        print("[RUNNING INFO]: Filtering GTF file\n")

    filtered_gtf = None

    if args.filterNonCoding:
        filtered_gtf = pbt.BedTool(my_gtf.filter(filterGtf, coding))
        feature_filtered = filtered_gtf.filter(
            filterFeatureGtf, ["exon"]
        ).saveas("filtered_gtf.gtf")
        filtered_gtf = pbt.BedTool(feature_filtered)
        print(pbt.BedTool(filtered_gtf).head())

    if args.filterCoding:
        if not args.featureTypes:
            sys.stderr.write(
                "Error: --filterCoding requires --featureTypes.\n"
            )
            sys.exit(-1)

        valid_features = {"exon", "gene", "transcript", "UTR", "CDS"}
        if not set(args.featureTypes).issubset(valid_features):
            sys.stderr.write(
                "Error: wrong featureType. Not in [exon, gene, transcript, UTR, CDS]!\n"
            )
            sys.exit(-1)

        filtered_gtf = pbt.BedTool(my_gtf.filter(filterGtf, non_coding))
        feature_filtered = filtered_gtf.filter(
            filterFeatureGtf, args.featureTypes
        ).saveas("filtered_gtf.gtf")
        filtered_gtf = pbt.BedTool(feature_filtered)
        print(pbt.BedTool(filtered_gtf).head())

    # ------------------------------------------------------------------ #
    # Intersect BED with filtered GTF                                     #
    # ------------------------------------------------------------------ #
    if filtered_gtf is not None:
        if args.verbose:
            print("[RUNNING INFO]: Running bedtools intersect on GTF and BED...\n")

        x = pbt.BedTool()
        intersect_bed = x.multi_intersect(
            i=[my_bed.fn, filtered_gtf.fn], header=True
        )
        print(pbt.BedTool(intersect_bed).head())

        df = pd.read_table(intersect_bed.fn, low_memory=False)
        df = df[df.num == 2]
        df.columns = df.iloc[0]
        df = df.iloc[:, [0, 1, 2]]
        intersect_path = args.oprefix + ".intersect.bed"
        df[1:].to_csv(intersect_path, sep="\t", index=False)

    # ------------------------------------------------------------------ #
    # Mosdepth (optional)                                                 #
    # ------------------------------------------------------------------ #
    if args.mosdepth:
        if args.verbose:
            print("[RUNNING INFO]: Running mosdepth...\n")

        sbp.run(
            [
                "mosdepth",
                "-t", str(args.thread),
                "--by", intersect_path,
                "-n",
                "--thresholds", str(args.minCoverage),
                "--mapq", str(args.minMapq),
                str(args.oprefix),
                str(args.bam),
            ]
        )
        sbp.run(["gunzip", args.oprefix + ".thresholds.bed.gz"])
        getEffGenomeSizeFromMosdepth(
            args.oprefix + ".thresholds.bed",
            use_mosdepth=True,
            verbose=True,
        )

        # Clean up mosdepth output files
        for suffix in [
            ".mosdepth.global.dist.txt",
            ".mosdepth.region.dist.txt",
            ".mosdepth.summary.txt",
            ".thresholds.bed.gz.csi",
            ".regions.bed.gz.csi",
        ]:
            _remove_if_exists(args.oprefix + suffix)

        if not args.saveIntermediates:
            _remove_if_exists(args.oprefix + ".regions.bed.gz")
            _remove_if_exists(args.oprefix + ".thresholds.bed")

    else:
        getEffGenomeSizeFromMosdepth(
            args.oprefix + ".intersect.bed",
            use_mosdepth=False,
            verbose=True,
        )

    print("\nScript Ended Successfully!\n")


def _remove_if_exists(path):
    """Remove *path* if it exists, silently ignore if not found."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
