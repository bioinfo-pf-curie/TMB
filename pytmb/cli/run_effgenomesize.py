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
- A genome annotation GTF file (.gtf or .gtf.gz)
- (optionally) a BAM file via mosdepth for coverage-based filtering

Typical usage::

    pyEffGenomeSize --bed design.bed --gtf genome.gtf.gz \\
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
# Biotype attribute strategy
# --------------------------
# GENCODE GTFs use:
#   gene-level entries    → gene_type
#   sub-gene entries      → transcript_type
# Ensembl GTFs use:
#   gene-level entries    → gene_biotype
#   sub-gene entries      → transcript_biotype
#
# filterGtf() checks the correct attribute depending on the feature column
# (col 3) so that both GENCODE and Ensembl files are handled correctly.
# ---------------------------------------------------------------------------

def filterGtf(interval, featuretype):
    """
    pybedtools filter callback: keep intervals whose biotype attribute is in
    *featuretype*.

    The biotype is resolved by checking attributes in the following priority
    order to support both GENCODE and Ensembl GTF conventions:

    * For **gene**-level entries (GTF feature column == ``gene``):
      ``gene_type`` (GENCODE) or ``gene_biotype`` (Ensembl).
    * For all other entries (``transcript``, ``exon``, ``CDS``, ``UTR`` …):
      ``transcript_type`` (GENCODE) or ``transcript_biotype`` (Ensembl).

    Parameters
    ----------
    interval : pybedtools.Interval
        A single GTF interval as parsed by pybedtools.
    featuretype : list of str
        Biotype values to keep (e.g. ``["protein_coding"]``).

    Returns
    -------
    str or False
        The matched biotype string when the interval passes the filter,
        ``False`` otherwise.
    """
    feature = interval[2]  # 3rd GTF column: gene / transcript / exon / CDS …

    if feature == "gene":
        # Gene-level entries use gene_type (GENCODE) or gene_biotype (Ensembl)
        biotype = (
            interval.attrs.get("gene_type")
            or interval.attrs.get("gene_biotype")
        )
    else:
        # Transcript / exon / CDS / UTR entries use transcript_type or transcript_biotype
        biotype = (
            interval.attrs.get("transcript_type")
            or interval.attrs.get("transcript_biotype")
        )

    if biotype and biotype in featuretype:
        return biotype
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

    # Required input
    parser.add_argument(
        "--bed", help="BED file (.bed)", required=True, default=None
    )

    # GTF options
    parser.add_argument(
        "--gtf", help="GTF file for genome annotation (.gtf or .gtf.gz)", default=None
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
            "(3rd GTF column).  Required with --filterCoding. Default: exon"
        ),
        nargs="+",
        default=["exon"],
    )

    # BAM options
    parser.add_argument(
        "--bam", help="BAM file for mapping statistics (.bam)", default=None
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

    # Output / misc
    parser.add_argument(
        "--saveIntermediates",
        help="Keep intermediate files (mosdepth output, filtered GTF, intersect BED)",
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
    # Mosdepth is automatically enabled when --bam is provided
    if args.bam:
        args.mosdepth = True
    elif args.minCoverage or args.minMapq:
        sys.stderr.write(
            "Error: --bam is required if --minCoverage or --minMapq is used.\n"
        )
        sys.exit(-1)

    if not (args.filterNonCoding or args.filterCoding):
        if args.gtf:
            sys.stderr.write(
                "Error: --filterCoding or --filterNonCoding is required when using GTF!\n"
            )
            sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Load BED / GTF                                                       #
    # ------------------------------------------------------------------ #
    if args.verbose:
        print("[INFO] Loading data...")

    my_bed = pbt.BedTool(args.bed)
    my_gtf = None
    filtered_gtf = None
    if args.gtf:
        my_gtf = pbt.BedTool(args.gtf)
        # Check that the GTF contains at least 3 entries and carries the
        # expected biotype attribute.  Both GENCODE (gene_type / transcript_type)
        # and Ensembl (gene_biotype / transcript_biotype) conventions are accepted.
        _BIOTYPE_ATTRS = {
            "gene_type", "gene_biotype",
            "transcript_type", "transcript_biotype",
        }
        try:
            third_entry = str(my_gtf[2])
        except IndexError:
            sys.stderr.write(
                "Error: GTF file has fewer than 3 entries. "
                "Cannot check for biotype attribute.\n"
            )
            sys.exit(-1)

        if not any(attr in third_entry for attr in _BIOTYPE_ATTRS):
            sys.stderr.write(
                "Error: GTF file has no recognised biotype attribute "
                "(gene_type, gene_biotype, transcript_type, transcript_biotype). "
                "Cannot filter this file.\n"
            )
            sys.exit(-1)

        # Count and warn about entries that carry none of the biotype attributes
        def _has_biotype(interval):
            return any(attr in interval.attrs for attr in _BIOTYPE_ATTRS)

        count = sum(1 for interval in my_gtf if not _has_biotype(interval))
        if count > 0:
            print(
                f"Warning: {count} annotations skipped because no biotype "
                "attribute (gene_type / transcript_type / *_biotype) was found."
            )

    # ------------------------------------------------------------------ #
    # Define coding / non-coding biotype sets                             #
    # Values cover both GENCODE and Ensembl naming conventions.           #
    # Legacy names (pre-GENCODE v24) are kept alongside current ones for  #
    # backwards compatibility.                                             #
    # ------------------------------------------------------------------ #
    coding = [
        # Protein-coding
        "protein_coding",
        # Immunoglobulin / T-cell receptor genes
        "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
        "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
        # Transcripts whose mRNA is degraded but originate from protein-coding loci
        "nonsense_mediated_decay", "non_stop_decay",
    ]
    non_coding = [
        # Long non-coding RNAs
        "lncRNA",            # current GENCODE label (was "lincRNA" before v24)
        "lincRNA",           # legacy label – retained for older GTFs
        "antisense_RNA",     # non-coding antisense transcript
        "antisense",         # legacy GENCODE label for antisense_RNA (pre-v24)
        "sense_intronic",
        "sense_overlapping",
        "3prime_overlapping_ncrna",
        "processed_transcript",
        "retained_intron",
        # Small / structural RNAs
        "miRNA", "misc_RNA",
        "Mt_rRNA", "Mt_tRNA",
        "rRNA", "ribozyme",
        "snoRNA", "snRNA",
        "scRNA", "scaRNA",
        "vault_RNA", "sRNA",
        # Pseudogenes
        "pseudogene",
        "polymorphic_pseudogene",
        "processed_pseudogene",
        "unprocessed_pseudogene",
        "unitary_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "transcribed_unitary_pseudogene",
        "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
        "IG_pseudogene",
        "TR_J_pseudogene", "TR_V_pseudogene",
    ]

    if args.verbose:
        print("[INFO] Filtering GTF file...")

    if my_gtf is not None:
        if args.filterNonCoding:
            filtered_gtf = pbt.BedTool(my_gtf.filter(filterGtf, coding))
            feature_filtered = filtered_gtf.filter(
                filterFeatureGtf, ["exon"]
            ).saveas("filtered_gtf.gtf")
            filtered_gtf = pbt.BedTool(feature_filtered)

        if args.filterCoding:
            if my_gtf is None:
                sys.stderr.write(
                    "Error: --filterCoding requires --gtf.\n"
                )
                sys.exit(-1)

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

    # ------------------------------------------------------------------ #
    # Intersect BED with filtered GTF                                     #
    # ------------------------------------------------------------------ #
    if filtered_gtf is not None:
        if args.verbose:
            print("[INFO] Running bedtools intersect...")

        # bedtools multiinter requires all inputs to be sorted by chrom then
        # start.  Sort both the capture BED and the filtered GTF before
        # calling multi_intersect.
        sorted_bed = my_bed.sort()
        sorted_gtf = filtered_gtf.sort()

        x = pbt.BedTool()
        intersect_bed = x.multi_intersect(
            i=[sorted_bed.fn, sorted_gtf.fn], header=True
        )

        df = pd.read_table(intersect_bed.fn, low_memory=False)
        df = df[df.num == 2]
        df.columns = df.iloc[0]
        df = df.iloc[:, [0, 1, 2]]
        intersect_path = args.oprefix + ".intersect.bed"
        df[1:].to_csv(intersect_path, sep="\t", index=False)

    # ------------------------------------------------------------------ #
    # Mosdepth (when BAM is provided)                                       #
    # ------------------------------------------------------------------ #
    if args.bam:
        if args.verbose:
            print("[INFO] Running mosdepth...")

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
        if not args.saveIntermediates:
            for suffix in [
                ".mosdepth.global.dist.txt",
                ".mosdepth.region.dist.txt",
                ".mosdepth.summary.txt",
                ".thresholds.bed.gz.csi",
                ".regions.bed.gz.csi",
            ]:
                _remove_if_exists(args.oprefix + suffix)

            _remove_if_exists(args.oprefix + ".regions.bed.gz")
            _remove_if_exists(args.oprefix + ".thresholds.bed")

    else:
        # Use BED file directly when GTF is not provided
        if filtered_gtf is None:
            getEffGenomeSizeFromMosdepth(
                args.bed,
                use_mosdepth=False,
                verbose=True,
            )
        else:
            getEffGenomeSizeFromMosdepth(
                args.oprefix + ".intersect.bed",
                use_mosdepth=False,
                verbose=True,
            )

    # ------------------------------------------------------------------ #
    # Cleanup temporary files                                              #
    # ------------------------------------------------------------------ #
    if not args.saveIntermediates:
        # Clean filtered GTF file
        _remove_if_exists("filtered_gtf.gtf")
        # Clean intersect BED file
        _remove_if_exists(args.oprefix + ".intersect.bed")

def _remove_if_exists(path):
    """Remove *path* if it exists, silently ignore if not found."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
