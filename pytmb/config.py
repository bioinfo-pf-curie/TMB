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
Configuration loading utilities for pyTMB.
"""

import yaml


def loadConfig(infile):
    """
    Load a YAML configuration file.

    Parameters
    ----------
    infile : str
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Parsed YAML content as a Python dictionary.
    """
    with open(infile, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError:
            raise
