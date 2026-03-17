#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of pyTMB software.
#
#  Copyright (c) 2020 - Institut Curie
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################
#
# COMPATIBILITY SHIM
# ------------------
# This script is kept for backward compatibility only.
# After installing the package with `pip install .` (or `pip install -e .`),
# prefer calling the `pyEffGenomeSize` command directly.
#
# All logic now lives in pytmb/cli/run_effgenomesize.py.
#
##############################################################################

from pytmb.cli.run_effgenomesize import main

if __name__ == "__main__":
    main()
