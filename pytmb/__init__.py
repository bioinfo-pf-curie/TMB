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

__version__ = "1.6.0"

from .config import loadConfig
from .vcf_utils import getTag, getMultiAlleleHeader
from .genome_size import getEffGenomeSizeFromBed
from .filters import (
    isAnnotatedAs,
    isPolym,
    isCancerHotspot,
    subsetINFO,
    infoTag2dl,
    info2dl,
)
from .tmb import calculate_tmb

__all__ = [
    "__version__",
    "loadConfig",
    "getTag",
    "getMultiAlleleHeader",
    "getEffGenomeSizeFromBed",
    "isAnnotatedAs",
    "isPolym",
    "isCancerHotspot",
    "subsetINFO",
    "infoTag2dl",
    "info2dl",
    "calculate_tmb",
]
