"""
PHACS-NMR - PHA Polymer Characterization by NMR.
Unified library for quantitative PHA copolymer composition, microstructural
statistics, and GUI analysis.

Copyright (C) 2026 Bioeconomy Science Institute
Licensed under the GNU General Public License v3.0 or later (GPLv3+).
See LICENSE for details.
"""

import sys

import numpy.linalg as _linalg

# Compatibility shim for legacy packages (e.g. nmrglue leastsqbound) expecting numpy.dual
sys.modules.setdefault("numpy.dual", _linalg)

from phacs_nmr import cli, core, gui
from phacs_nmr.cli import (
    characteristic_peaks,
    find_polymer_ID,
    main,
    parse_args,
    start_analysis,
)
from phacs_nmr.core.io import treatment_of_txt_file
from phacs_nmr.core.noise import calc_noise_level, inv_scale
from phacs_nmr.path_utils import (
    is_drive_or_root,
    resolve_output_prefix,
    resolve_sample_name,
)

__version__ = "2.0.0"

__all__ = [
    "__version__",
    "calc_noise_level",
    "characteristic_peaks",
    "cli",
    "core",
    "find_polymer_ID",
    "gui",
    "inv_scale",
    "is_drive_or_root",
    "main",
    "parse_args",
    "resolve_output_prefix",
    "resolve_sample_name",
    "start_analysis",
    "treatment_of_txt_file",
]
