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
Core TMB calculation logic for pyTMB.

The public entry-point is :func:`calculate_tmb`, which accepts all filtering
parameters and returns a result dict.  The CLI wrappers in
:mod:`pytmb.cli.run_tmb` call this function and handle I/O.
"""

import re
import sys
import warnings
import os.path
from datetime import date

import cyvcf2

from .vcf_utils import getTag
from .filters import (
    isAnnotatedAs,
    isPolym,
    isCancerHotspot,
    infoTag2dl,
    info2dl,
)


def calculate_tmb(
    vcf_path,
    db_flags,
    caller_flags,
    eff_genome_size,
    sample=None,
    vaf=0,
    maf=1,
    min_depth=1,
    min_alt_depth=1,
    filter_low_qual=False,
    filter_indels=False,
    filter_coding=False,
    filter_splice=False,
    filter_non_coding=False,
    filter_syn=False,
    filter_non_syn=False,
    filter_cancer_hotspot=False,
    filter_polym=False,
    filter_recurrence=False,
    polym_db="gnomad",
    cancer_db="cosmic",
    export_path=None,
    debug=False,
    verbose=False,
):
    """
    Run the TMB calculation on a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the input VCF/BCF file.
    db_flags : dict
        Parsed database configuration (from a ``dbConfig`` YAML file).
    caller_flags : dict
        Parsed variant-caller configuration (from a ``varConfig`` YAML file).
    eff_genome_size : int
        Effective genome size in base pairs.
    sample : str or None, optional
        Sample ID to focus on.  Required when the VCF contains more than one
        sample.
    vaf : float, optional
        Filter variants with allelic ratio < *vaf*. Default 0.
    maf : float, optional
        Filter variants with MAF > *maf*. Default 1.
    min_depth : int, optional
        Filter variants with total depth < *min_depth*. Default 1.
    min_alt_depth : int, optional
        Filter variants with alt-allele depth < *min_alt_depth*. Default 1.
    filter_low_qual : bool, optional
        Discard variants that do not have PASS in the FILTER column.
    filter_indels : bool, optional
        Discard insertion/deletion variants.
    filter_coding : bool, optional
        Discard coding variants.
    filter_splice : bool, optional
        Discard splice-site variants.
    filter_non_coding : bool, optional
        Discard non-coding variants.
    filter_syn : bool, optional
        Discard synonymous variants.
    filter_non_syn : bool, optional
        Discard non-synonymous variants.
    filter_cancer_hotspot : bool, optional
        Discard variants annotated as cancer hotspots.
    filter_polym : bool, optional
        Discard variants present in population databases above *maf*.
    filter_recurrence : bool, optional
        Discard variants flagged by run-level recurrence.
    polym_db : str, optional
        Comma-separated list of population databases to use (e.g.
        ``"gnomad,1k"``).  Default ``"gnomad"``.
    cancer_db : str, optional
        Comma-separated list of cancer hotspot databases (e.g. ``"cosmic"``).
        Default ``"cosmic"``.
    export_path : str or None, optional
        If provided, write the variants that pass all filters to this path as
        a compressed VCF (``wz`` mode).
    debug : bool, optional
        Write every variant to a ``_debug.vcf.gz`` file with an added
        ``TMB_FILTERS`` INFO tag listing which filters were triggered.
    verbose : bool, optional
        Print progress messages every 1 000 variants.

    Returns
    -------
    dict
        A result dict with the following keys:

        * ``tmb`` – float, computed TMB (mutations / Mb)
        * ``var_counter`` – int, total variants seen
        * ``var_ni`` – int, variants skipped due to missing annotation
        * ``var_tmb`` – int, variants counted toward TMB
        * ``eff_genome_size`` – int, effective genome size used
        * ``filter_stats`` – dict, per-filter counts

    Raises
    ------
    SystemExit
        On critical input errors (sample not found, empty genome size, etc.).
    """
    # ------------------------------------------------------------------ #
    # Open VCF                                                             #
    # ------------------------------------------------------------------ #
    if sample is not None:
        vcf = cyvcf2.VCF(vcf_path, samples=sample)
        if sample not in vcf.samples:
            sys.stderr.write(f"Error: Sample '{sample}' not found in VCF\n")
            sys.exit(-1)
    else:
        vcf = cyvcf2.VCF(vcf_path)

    if len(vcf.samples) > 1:
        sys.stderr.write(
            "Error: "
            + str(len(vcf.samples))
            + " samples detected. This version is designed for a single sample!"
            " Use --sample argument.\n"
        )
        sys.exit(-1)

    # ------------------------------------------------------------------ #
    # Optional output writers                                              #
    # ------------------------------------------------------------------ #
    wx = None
    wd = None

    if export_path is not None:
        wx = cyvcf2.Writer(export_path, vcf, mode="wz")

    if debug:
        vcf.add_info_to_header(
            {
                "ID": "TMB_FILTERS",
                "Description": "Detected filters for TMB calculation",
                "Type": "Character",
                "Number": "1",
            }
        )
        debug_path = re.sub(
            r"\.vcf$|\.vcf\.gz$|\.bcf", "_debug.vcf.gz", os.path.basename(vcf_path)
        )
        wd = cyvcf2.Writer(debug_path, vcf, mode="wz")

    # ------------------------------------------------------------------ #
    # Pre-compute database field lists (avoid rebuilding per-variant)      #
    # ------------------------------------------------------------------ #
    cancer_fields = []
    if filter_cancer_hotspot:
        for db in cancer_db.split(","):
            cancer_fields.extend(db_flags["cancerDb"][db])

    polym_fields = []
    if filter_polym:
        for db in polym_db.split(","):
            polym_fields.extend(db_flags["polymDb"][db])

    # ------------------------------------------------------------------ #
    # Counters                                                             #
    # ------------------------------------------------------------------ #
    var_counter = 0
    var_ni = 0
    var_tmb = 0

    filter_stats = {
        "INDEL": 0,
        "QUAL": 0,
        "VAF": 0,
        "DEPTH": 0,
        "ALTDEPTH": 0,
        "CODING": 0,
        "SPLICING": 0,
        "NONCODING": 0,
        "SYN": 0,
        "NON_SYN": 0,
        "HOTSPOT": 0,
        "POLYM": 0,
        "RUNREC": 0,
    }

    # ------------------------------------------------------------------ #
    # Main variant loop                                                    #
    # ------------------------------------------------------------------ #
    for variant in vcf:
        var_counter += 1
        if verbose and var_counter % 1000 == 0:
            print("## ", var_counter)
            if debug and var_counter == 1000:
                sys.exit()

        try:
            db_info = dict(variant.INFO)
            debug_info = ""

            # Parse annotation INFO tag
            if db_flags["tag"] != "":
                annot_info = infoTag2dl(variant.INFO.get(db_flags["tag"]))
            else:
                annot_info = info2dl(variant.INFO)

            # Skip variants with no annotations
            if db_info is None or annot_info is None:
                var_ni += 1
                continue

            # ---- Indels ------------------------------------------------
            if filter_indels and variant.is_indel:
                debug_info = ",".join([debug_info, "INDEL"])
                filter_stats["INDEL"] += 1
                if not debug:
                    continue

            # ---- Low quality -------------------------------------------
            if filter_low_qual and variant.FILTER is not None:
                debug_info = ",".join([debug_info, "QUAL"])
                filter_stats["QUAL"] += 1
                if not debug:
                    continue

            # ---- Variant Allele Frequency --------------------------------
            fval = getTag(variant, caller_flags["freq"])
            if fval is not None and len(fval[fval < vaf]) == len(variant.ALT):
                debug_info = ",".join([debug_info, "VAF"])
                filter_stats["VAF"] += 1
                if not debug:
                    continue

            # ---- Sequencing Depth ----------------------------------------
            dval = getTag(variant, caller_flags["depth"])
            if dval is not None and len(dval[dval < min_depth]) == len(variant.ALT):
                debug_info = ",".join([debug_info, "DEPTH"])
                filter_stats["DEPTH"] += 1
                if not debug:
                    continue

            # ---- Alternative Allele Depth --------------------------------
            ad = getTag(variant, caller_flags["altDepth"]).flatten()
            # If AD = REF + ALTs, remove the first (REF) value
            if len(ad) == (len(variant.ALT) + 1):
                ad = ad[1:]
            if ad is not None and len(ad[ad < min_alt_depth]) == len(variant.ALT):
                debug_info = ",".join([debug_info, "ALTDEPTH"])
                filter_stats["ALTDEPTH"] += 1
                if not debug:
                    continue

            # ---- Coding variants -----------------------------------------
            if filter_coding and isAnnotatedAs(
                variant, infos=annot_info, flags=db_flags["isCoding"], sep=db_flags["sep"]
            ):
                if not isAnnotatedAs(
                    variant,
                    infos=annot_info,
                    flags=db_flags["isNonCoding"],
                    sep=db_flags["sep"],
                ):
                    debug_info = ",".join([debug_info, "CODING"])
                    filter_stats["CODING"] += 1
                    if not debug:
                        continue

            # ---- Splice variants -----------------------------------------
            if filter_splice and isAnnotatedAs(
                variant,
                infos=annot_info,
                flags=db_flags["isSplicing"],
                sep=db_flags["sep"],
            ):
                debug_info = ",".join([debug_info, "SPLICING"])
                filter_stats["SPLICING"] += 1
                if not debug:
                    continue

            # ---- Non-coding variants -------------------------------------
            if filter_non_coding and isAnnotatedAs(
                variant,
                infos=annot_info,
                flags=db_flags["isNonCoding"],
                sep=db_flags["sep"],
            ):
                if not isAnnotatedAs(
                    variant,
                    infos=annot_info,
                    flags=db_flags["isCoding"],
                    sep=db_flags["sep"],
                ):
                    debug_info = ",".join([debug_info, "NONCODING"])
                    filter_stats["NONCODING"] += 1
                    if not debug:
                        continue

            # ---- Synonymous variants -------------------------------------
            if filter_syn and isAnnotatedAs(
                variant,
                infos=annot_info,
                flags=db_flags["isSynonymous"],
                sep=db_flags["sep"],
            ):
                if not isAnnotatedAs(
                    variant,
                    infos=annot_info,
                    flags=db_flags["isNonSynonymous"],
                    sep=db_flags["sep"],
                ):
                    debug_info = ",".join([debug_info, "SYN"])
                    filter_stats["SYN"] += 1
                    if not debug:
                        continue

            # ---- Non-synonymous variants ---------------------------------
            if filter_non_syn and isAnnotatedAs(
                variant,
                infos=annot_info,
                flags=db_flags["isNonSynonymous"],
                sep=db_flags["sep"],
            ):
                if not isAnnotatedAs(
                    variant,
                    infos=annot_info,
                    flags=db_flags["isSynonymous"],
                    sep=db_flags["sep"],
                ):
                    debug_info = ",".join([debug_info, "NON_SYN"])
                    filter_stats["NON_SYN"] += 1
                    if not debug:
                        continue

            # ---- Cancer hotspot ------------------------------------------
            if filter_cancer_hotspot:
                if isCancerHotspot(variant, infos=db_info, flags=cancer_fields):
                    debug_info = ",".join([debug_info, "HOTSPOT"])
                    filter_stats["HOTSPOT"] += 1
                    if not debug:
                        continue

            # ---- Polymorphisms -------------------------------------------
            if filter_polym:
                if isPolym(variant, infos=db_info, flags=polym_fields, val=maf):
                    debug_info = ",".join([debug_info, "POLYM"])
                    filter_stats["POLYM"] += 1
                    if not debug:
                        continue

            # ---- Run-level recurrence ------------------------------------
            if filter_recurrence:
                if isPolym(
                    variant, infos=db_info, flags=db_flags["recurrence"]["run"], val=0
                ):
                    debug_info = ",".join([debug_info, "RUNREC"])
                    filter_stats["RUNREC"] += 1
                    if not debug:
                        continue

        except Exception:
            warn_flag = (
                str(variant.CHROM) + ":" + str(variant.start) + "-" + str(variant.end)
            )
            warnings.warn(
                f"Warning: variant {warn_flag} raises an error. Skipping."
            )
            continue

        # ------------------------------------------------------------------
        # Variant passes all active filters
        # ------------------------------------------------------------------
        if debug_info == "":
            var_tmb += 1
            if wx is not None:
                wx.write_record(variant)

        if debug and wd is not None:
            variant.INFO["TMB_FILTERS"] = re.sub(r"^,", "", debug_info)
            wd.write_record(variant)

    # ------------------------------------------------------------------ #
    # Close writers                                                        #
    # ------------------------------------------------------------------ #
    if wx is not None:
        wx.close()
    if wd is not None:
        wd.close()
    vcf.close()

    # ------------------------------------------------------------------ #
    # Compute TMB                                                          #
    # ------------------------------------------------------------------ #
    tmb = round(float(var_tmb) / (float(eff_genome_size) / 1e6), 2)

    return {
        "tmb": tmb,
        "var_counter": var_counter,
        "var_ni": var_ni,
        "var_tmb": var_tmb,
        "eff_genome_size": eff_genome_size,
        "filter_stats": filter_stats,
    }


def print_tmb_report(results, args, version):
    """
    Print the TMB summary report to stdout.

    Parameters
    ----------
    results : dict
        The dict returned by :func:`calculate_tmb`.
    args : argparse.Namespace
        Parsed CLI arguments (used to echo back the filter settings).
    version : str
        Package version string.
    """
    var_counter = results["var_counter"]
    filter_stats = results["filter_stats"]

    if results["var_ni"] == var_counter:
        print("Warning: No Annotation detected. Is your file annotated ?")
        sys.exit(-1)

    print("pyTMB version=", version)
    print("When=", date.today())
    print("")
    print("Input=", args.vcf)
    print("Sample=", args.sample)
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
    print("filterCancerHotspot=", args.filterCancerHotspot)
    print("filterPolym=", args.filterPolym)
    print("")
    print("Filter statistics:")
    print("------------------")
    for filt, count in filter_stats.items():
        if count > 0:
            pct = round(100.0 * count / var_counter, 2) if var_counter > 0 else 0
            print(f"  {filt}= {count} ({pct}%)")
    print("")
    print("Total number of variants=", var_counter)
    print("Non-informative variants=", results["var_ni"])
    print("Variants after filters=", results["var_tmb"])
    print("Effective Genome Size=", results["eff_genome_size"])
    print("")
    print("TMB=", results["tmb"])
