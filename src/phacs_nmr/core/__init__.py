"""
PHACS-NMR core engine package.
Provides modular NMR processing, statistical copolymer microstructure models,
declarative polymer profiles, and report generators.

Copyright (C) 2026 Bioeconomy Science Institute
Licensed under the GNU General Public License v3.0 or later (GPLv3+).
See LICENSE for details.
"""

from phacs_nmr.core.config import (
    P34HB_PROFILE,
    PHBHHX_PROFILE,
    PHBHHXHO_MCL_PROFILE,
    PHBV_PROFILE,
    PROFILES,
    PolymerProfile,
    get_profile,
)
from phacs_nmr.core.engine import (
    compute_area_and_ratio,
    compute_area_and_ratio_monomers,
    compute_area_and_ratio_monomers_mcl,
    create_pdf_report,
    export_csv_results,
    run_analysis,
)
from phacs_nmr.core.io import (
    data_dcon,
    find_monomers,
    read_ascii_spectrum,
    read_dcon_file,
    results_dconv,
    treatment_of_txt_file,
)
from phacs_nmr.core.models import (
    apply_copolymer_statistics,
    calculate_average_lengths,
    calculate_chi_parameter,
    calculate_d_parameter,
    compute_kamiya_model_1,
    compute_kamiya_model_2,
    compute_kamiya_model_3,
)
from phacs_nmr.core.noise import calc_noise_level, inv_scale
from phacs_nmr.core.peak_fitting import split_graphic_regions

__all__ = [
    "P34HB_PROFILE",
    "PHBHHXHO_MCL_PROFILE",
    "PHBHHX_PROFILE",
    "PHBV_PROFILE",
    "PROFILES",
    "PolymerProfile",
    "apply_copolymer_statistics",
    "calc_noise_level",
    "calculate_average_lengths",
    "calculate_chi_parameter",
    "calculate_d_parameter",
    "compute_area_and_ratio",
    "compute_area_and_ratio_monomers",
    "compute_area_and_ratio_monomers_mcl",
    "compute_kamiya_model_1",
    "compute_kamiya_model_2",
    "compute_kamiya_model_3",
    "create_pdf_report",
    "data_dcon",
    "export_csv_results",
    "find_monomers",
    "get_profile",
    "inv_scale",
    "read_ascii_spectrum",
    "read_dcon_file",
    "results_dconv",
    "run_analysis",
    "split_graphic_regions",
    "treatment_of_txt_file",
]
