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
VCF utility functions for pyTMB.
"""

import numpy as np


def getMultiAlleleHeader(vcf):
    """
    Extract VCF column names that carry multi-allelic information
    (i.e. fields annotated with Number='A' in the VCF header).

    Parameters
    ----------
    vcf : cyvcf2.VCF
        An open cyvcf2 VCF reader.

    Returns
    -------
    dict
        A dict with keys 'FORMAT' and 'INFO', each containing a list of
        field IDs that have Number='A'.
    """
    FORMAT = []
    INFO = []
    for h in vcf.header_iter():
        i = h.info(extra=True)
        if "Number" in i.keys() and i["Number"] == "A":
            if i["HeaderType"] == "FORMAT":
                FORMAT.append(i["ID"])
            elif i["HeaderType"] == "INFO":
                INFO.append(i["ID"])
    return dict(FORMAT=FORMAT, INFO=INFO)


def getTag(v, tag):
    """
    Retrieve the value of *tag* from a cyvcf2 variant record.

    The FORMAT field is checked first; if the tag is absent or returns None
    there, the INFO field is tried.  The result is always returned as a
    1-D float ``numpy.ndarray``.

    Parameters
    ----------
    v : cyvcf2.Variant
        A variant record from a cyvcf2 VCF iterator.
    tag : str
        The FORMAT or INFO field name to retrieve.

    Returns
    -------
    numpy.ndarray
        A 1-D float array with the tag value(s).
    """
    # First check in FORMAT field
    val = None
    if tag in v.FORMAT:
        val = v.format(tag)

    # Otherwise, check in INFO field
    if tag not in v.FORMAT or val is None:
        val = v.INFO.get(tag)

    if not isinstance(val, np.ndarray):
        val = np.array([val], float)
    else:
        val = val.astype("float")

    return val
