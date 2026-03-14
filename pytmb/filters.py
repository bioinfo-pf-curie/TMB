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
Variant annotation and filtering functions for pyTMB.
"""


def subsetINFO(annot, keys):
    """
    Subset an annotation dict (or list of dicts) to only the provided keys.

    Parameters
    ----------
    annot : dict or list of dict
        The full annotation information for a variant (e.g. from INFO).
    keys : iterable
        The keys to keep.

    Returns
    -------
    dict or list of dict
        Filtered annotation containing only the requested keys.
    """
    if isinstance(annot, list):
        subsetInfo = []
        for i in range(len(annot)):
            z = {k: annot[i][k] for k in keys if k in annot[i]}
            if len(z) > 0:
                subsetInfo.append(z)
    else:
        subsetInfo = {k: annot[k] for k in keys if k in annot}
    return subsetInfo


def infoTag2dl(INFO):
    """
    Parse a pipe-delimited annotation tag value (e.g. snpEff ANN field)
    into a list of dicts keyed by positional index.

    Parameters
    ----------
    INFO : str or None
        Raw INFO tag value (comma-separated annotations, each pipe-separated).

    Returns
    -------
    list of dict or None
        List of annotation dicts, or None if *INFO* is None.
    """
    if INFO is not None:
        annotTag = INFO.split(",")
        annotInfo = []
        for tag in annotTag:
            annot = tag.split("|")
            dictannot = {i: annot[i] for i in range(len(annot))}
            annotInfo.append(dictannot)
        return annotInfo
    return None


def info2dl(INFO):
    """
    Wrap an INFO dict into a one-element list of dicts (for ANNOVAR-style
    flat INFO fields).

    Parameters
    ----------
    INFO : mapping or None
        A cyvcf2 INFO mapping.

    Returns
    -------
    list of dict or None
        ``[dict(INFO)]``, or None if *INFO* is None.
    """
    if INFO is not None:
        return [dict(INFO)]
    return None


def isAnnotatedAs(v, infos, flags, sep):
    """
    Check whether a variant carries any of the given annotation flags.

    Parameters
    ----------
    v : cyvcf2.Variant
        The variant record (not used directly, kept for API consistency).
    infos : list of dict
        Parsed annotation information (output of :func:`infoTag2dl` or
        :func:`info2dl`).
    flags : dict
        Mapping of annotation field key → list of expected values.
    sep : str
        Separator character used between multiple effect annotations
        within a single field value.

    Returns
    -------
    bool
        True if any annotation in *infos* matches any value in *flags*.
    """
    subINFO = subsetINFO(infos, keys=flags.keys())
    for sub in subINFO:
        for k in flags.keys():
            for val in flags[k]:
                for subval in sub[k].split(sep):
                    if val == subval:
                        return True
    return False


def isPolym(v, infos, flags, val):
    """
    Check whether a variant is present in a population genome database with
    a minor-allele frequency (MAF) at or above *val*.

    Parameters
    ----------
    v : cyvcf2.Variant
        The variant record (not used directly, kept for API consistency).
    infos : dict
        Flat INFO dict for the variant (``dict(variant.INFO)``).
    flags : list of str
        INFO field names that carry population-frequency values.
    val : float
        MAF threshold; variants at or above this frequency are considered
        polymorphisms.

    Returns
    -------
    bool
        True if any population frequency field meets or exceeds *val*.
    """
    subINFO = subsetINFO(infos, keys=flags)
    for key in subINFO:
        if isinstance(subINFO[key], tuple):
            for i in subINFO[key]:
                if i is not None and i != ".":
                    if float(i) >= float(val):
                        return True
        elif (
            subINFO[key] is not None
            and subINFO[key] != "."
            and subINFO[key] != "NA"
        ):
            if float(subINFO[key]) >= float(val):
                return True
    return False


def isCancerHotspot(v, infos, flags):
    """
    Check whether a variant is annotated as a cancer hotspot.

    Parameters
    ----------
    v : cyvcf2.Variant
        The variant record (not used directly, kept for API consistency).
    infos : dict
        Flat INFO dict for the variant (``dict(variant.INFO)``).
    flags : list of str
        INFO field names used to flag cancer hotspots.

    Returns
    -------
    bool
        True if any hotspot INFO field is present and non-empty.
    """
    subINFO = subsetINFO(infos, keys=flags)
    for key in subINFO:
        if subINFO[key] is not None and subINFO[key] != ".":
            return True
    return False
