"""
PHACS-NMR - PHA Polymer Characterization by NMR.
Unified entry point and command-line execution runner.

Option A Smart Launcher:
- No arguments or --gui: launches the Graphical User Interface (GUI).
- File path provided: runs single-spectrum polymer auto-detection and deconvolution.
- Directory path provided: runs batch headless processing on all *.txt spectra.

Copyright (C) 2026 Bioeconomy Science Institute
Licensed under the GNU General Public License v3.0 or later (GPLv3+).
See LICENSE for details.
"""
import argparse
import sys
from pathlib import Path
from typing import Any

import nmrglue as ng
import numpy as np

from phacs_nmr.core.config import (
    P34HB_PROFILE,
    PHBHHX_PROFILE,
    PHBHHXHO_MCL_PROFILE,
    PHBV_PROFILE,
)
from phacs_nmr.core.engine import run_analysis
from phacs_nmr.core.io import treatment_of_txt_file
from phacs_nmr.core.noise import calc_noise_level, inv_scale


def characteristic_peaks(ppm_scale: Any, data: Any) -> list[list[Any]]:
    """Find specific diagnostic peaks for each monomer on the spectrum."""
    limits = []
    limits.append(['X6', inv_scale(ppm_scale, 13.7), inv_scale(ppm_scale, 13.9)])
    limits.append(['B4', inv_scale(ppm_scale, 19.5), inv_scale(ppm_scale, 20)])
    limits.append(['O7', inv_scale(ppm_scale, 22.4), inv_scale(ppm_scale, 22.7)])
    limits.append(['F2', inv_scale(ppm_scale, 30.5), inv_scale(ppm_scale, 30.9)])
    limits.append(['V5', inv_scale(ppm_scale, 9.2), inv_scale(ppm_scale, 9.5)])
    limits.append(['mcl', inv_scale(ppm_scale, 25), inv_scale(ppm_scale, 25.2)])
    return limits


def find_polymer_ID(ppm_scale: Any, data: Any, multiple_noise: int) -> list[int]:
    """Detect polymer composition ID indicating the constitutive monomers present."""
    noise = calc_noise_level(ppm_scale, data)
    threshold = int(multiple_noise * noise)
    lineshapes = ['l']
    peak_table = ng.peakpick.pick(
        np.array(data),
        pthres=threshold,
        algorithm='downward',
        lineshapes=lineshapes,
    )
    peaks = peak_table
    limits = characteristic_peaks(ppm_scale, data)
    polymer_id = [0, 0, 0, 0, 0, 0]

    for _id, info in enumerate(limits):
        start, end = min(info[1], info[2]), max(info[1], info[2])
        params = [[(x, max(0.1, lw))] for x, lw in zip(peaks['X_AXIS'], peaks['X_LW']) if start < x < end]
        npeaks = len(params)
        if info[0] == 'F2' and npeaks != 0:
            polymer_id[0] = 1
        if info[0] == 'B4' and npeaks != 0:
            polymer_id[1] = 1
        if info[0] == 'V5' and npeaks != 0:
            polymer_id[2] = 1
        if info[0] == 'X6' and npeaks != 0:
            polymer_id[3] = 1
        if info[0] == 'O7' and npeaks != 0:
            polymer_id[4] = 1
        if info[0] == 'mcl' and npeaks != 0:
            polymer_id[5] = 1

    return polymer_id


def start_analysis(input_file: str, multiple_noise: int = 10) -> list[Any]:
    """
    Run automated polymer identification, peak deconvolution, and reporting
    for a single 13C ASCII spectrum file.
    """
    print('Working on: ' + input_file)
    print('Treatment of NMR txt file')
    x, y = treatment_of_txt_file(input_file)
    print('text file treated')
    p_ID = find_polymer_ID(x, y, multiple_noise)
    code = ['F', 'B', 'V', 'X', 'O', 'mcl']

    if p_ID == [0, 1, 0, 1, 1, 0]:
        print("Polymer type: PHBHHxHO")
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 0, 0, 1, 1, 0]:
        print("Polymer type: PHHxHO")
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 0, 0, 1, 0]:
        print("Polymer type: PHBHO", p_ID)
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 0, 0, 1, 1]:
        print("Polymer type: PHBHO_mcl", p_ID)
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 0, 1, 1, 1]:
        print("Polymer type: PHBHHxHO_mcl")
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 0, 0, 1, 1, 1]:
        print("Polymer type: PHHxHO_mcl")
        figures = run_analysis(input_file, PHBHHXHO_MCL_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [1, 1, 0, 0, 0, 0]:
        print("Polymer type: P34HB")
        figures = run_analysis(input_file, P34HB_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 0, 1, 0, 0]:
        print("Polymer type: PHBHHx")
        print('Job done')
        figures = run_analysis(input_file, PHBHHX_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 1, 0, 0, 0]:
        print("Polymer type: PHBV")
        print("Job done")
        figures = run_analysis(input_file, PHBV_PROFILE, multiple_noise=multiple_noise)
    elif p_ID == [0, 1, 0, 0, 0, 0]:
        figures = []
        print("Polymer type: PHB homopolymer")
    else:
        figures = []
        print("Polymer ID:", p_ID, code)
        print("The model is not enough developed for other complex heteropolymers")

    return figures


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="PHACS-NMR - PHA Polymer Characterization by NMR (GUI & CLI)."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=None,
        help="Path to an input NMR ASCII .txt file (single spectrum) or directory (batch processing). If omitted, launches GUI.",
    )
    parser.add_argument(
        "noise_level",
        type=int,
        nargs="?",
        default=10,
        choices=(5, 10),
        help="Noise multiplier (5 for clean signal, 10 for noisy signal; default: 10).",
    )
    parser.add_argument(
        "-n",
        "--noise",
        type=int,
        dest="opt_noise",
        choices=(5, 10),
        default=None,
        help="Noise multiplier option (5 or 10; overrides positional noise_level if specified).",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Explicitly launch the Graphical User Interface.",
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> Any:
    """Primary execution entrypoint supporting GUI launch, single spectrum, or batch directory processing."""
    parsed = parse_args(args)

    if parsed.gui or parsed.path is None:
        from phacs_nmr.gui import main as gui_main
        return gui_main()

    target = Path(parsed.path).resolve()
    if not target.exists():
        raise FileNotFoundError(f"Input path not found: {target}")

    noise = parsed.opt_noise if parsed.opt_noise is not None else parsed.noise_level

    current_mod = sys.modules[__name__]
    if target.is_file():
        return current_mod.start_analysis(str(target), noise)
    elif target.is_dir():
        txt_files = sorted(target.glob("*.txt"))
        if not txt_files:
            print(f"No .txt files found in {target}")
            return []
        results = []
        for input_file in txt_files:
            print(f"\n=== Processing {input_file.name} ===")
            try:
                res = current_mod.start_analysis(str(input_file), noise)
                results.append((input_file, res))
            except Exception as exc:
                print(f"Failed on {input_file}: {exc}")
        return results
    else:
        raise FileNotFoundError(f"Unsupported path type: {target}")


if __name__ == "__main__":
    main()
