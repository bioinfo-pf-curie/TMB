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
Effective genome size computation utilities for pyTMB.
"""

import sys


def getEffGenomeSizeFromBed(infile, verbose=False):
    """
    Calculate the effective genome size (in bp) from a BED file by summing
    the lengths of all intervals.

    Parameters
    ----------
    infile : str
        Path to the BED file.
    verbose : bool, optional
        Print progress messages when True. Default is False.

    Returns
    -------
    int
        Total effective genome size in base pairs.

    Raises
    ------
    SystemExit
        On malformed input or missing file.
    """
    if verbose:
        print(f"## File Name: {infile}")

    effgs = 0
    nline = 0

    try:
        with open(infile) as bedhandle:
            for line in bedhandle:
                nline += 1
                if line.startswith("#") or line.strip() == "":
                    continue  # skip comments or empty lines

                bedtab = line.strip().split("\t")
                if len(bedtab) < 3:
                    sys.stderr.write(
                        f"Error: wrong input format in line {nline}. Not a BED file?\n"
                    )
                    sys.exit(-1)

                try:
                    start = int(bedtab[1])
                    end = int(bedtab[2])
                except ValueError:
                    sys.stderr.write(
                        f"Error: non-integer coordinates in line {nline}.\n"
                    )
                    sys.exit(-1)

                effgs += abs(end - start)
    except FileNotFoundError:
        sys.stderr.write(f"Error: file '{infile}' not found.\n")
        sys.exit(-1)

    if verbose:
        print(f"## Effective Genome Size = {effgs} bp")

    return effgs


def getEffGenomeSizeFromMosdepth(infile, use_mosdepth=False, verbose=False):
    """
    Calculate the effective genome size from a mosdepth thresholds BED file
    (or a plain BED file).

    When *use_mosdepth* is True the file is expected to have a header line and
    a 5th column containing the per-region coverage value; the callable region
    size is reported as the sum of those coverage values.

    Parameters
    ----------
    infile : str
        Path to the mosdepth thresholds BED (or plain BED) file.
    use_mosdepth : bool, optional
        Set to True when processing a mosdepth thresholds file.
        Default is False.
    verbose : bool, optional
        Print summary statistics when True. Default is False.

    Returns
    -------
    int
        Effective (callable) genome size.

    Raises
    ------
    SystemExit
        On malformed input.
    """
    effgs = 0
    totgs = 0
    nline = 0

    with open(infile, "rt") as f:
        if use_mosdepth:
            next(f)  # skip header

        for line in f:
            bedtab = line.strip().split("\t")
            nline += 1

            if use_mosdepth:
                try:
                    chromosome, start, end, region, coverage = bedtab[:5]
                    effgs += int(coverage)
                except ValueError:
                    sys.stderr.write(
                        f"Error: wrong input format in line {nline}. Not a BED file?\n"
                    )
                    sys.exit(-1)
            else:
                try:
                    chromosome, start, end = bedtab[:3]
                except ValueError:
                    sys.stderr.write(
                        f"Error: wrong input format in line {nline}. Not a BED file?\n"
                    )
                    sys.exit(-1)

            intl = abs(int(end) - int(start))
            totgs += intl

    if verbose:
        print(f"## File Name: {infile}")
        print(f"## Total region size = {totgs}")
        if use_mosdepth:
            pct = round(effgs / totgs * 100, 3) if totgs > 0 else 0
            print(f"## Callable region = {effgs} ({pct}%)")
        print(f"## Effective Genome Size = {totgs if not use_mosdepth else effgs} bp\n")

    return effgs if use_mosdepth else totgs
