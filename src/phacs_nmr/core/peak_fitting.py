"""
Spectral region segmentation and NMR lineshape fitting using nmrglue.
Deconstructs NMR spectra into diagnostic chemical shift windows and fits
pseudo-Voigt/Lorentzian peaks to determine resonance areas and peak heights.

References:
    Helmus, J.J.; Jaroniec, C.P. nmrglue: An open source Python package for
    the analysis of multidimensional NMR data. J. Biomol. NMR 2013, 55, 355-367.
    DOI: 10.1007/s10858-013-9718-x.
"""
from collections.abc import Callable
from math import ceil, pi, sqrt
from typing import Any

import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np
from matplotlib.figure import Figure

from phacs_nmr.core.io import read_ascii_spectrum
from phacs_nmr.core.noise import calc_noise_level, inv_scale


def split_graphic_regions(
    file_path: str | Any,
    table: list[tuple[float, float, str]],
    multiple_noise: float,
    limits_fn: Callable[[Any], tuple[list[list[Any]], dict[str, tuple[int, int]]]],
    col_x_index: int = 3,
    col_y_index: int = 1,
    col_x_name: str = "ppm",
    col_y_name: str = "Intensity",
    header: Any = None,
    figsize: tuple[float, float] = (10.0, 7.0),
    noise_as_int: bool = False,
    monomer_colors: dict[str, str] | None = None,
    preloaded_xy: tuple[list[float], list[float]] | None = None,
) -> tuple[list[float], list[float], dict[str, Any], Figure]:
    """
    Split NMR spectrum into narrower regions to optimize peak identification
    and create subplots for each region.
    """
    n_rows = ceil(len(table) / 4)
    fig, axes = plt.subplots(n_rows, 4, squeeze=False, figsize=figsize)

    x: list[float] = []
    y: list[float] = []
    pha_df: dict[str, Any] = {}

    try:
        if preloaded_xy is not None:
            x, y = preloaded_xy
        else:
            x, y = read_ascii_spectrum(file_path, col_x=col_x_index, col_y=col_y_index)

        raw_noise = multiple_noise * calc_noise_level(x, y)
        noise = int(raw_noise) if noise_as_int else raw_noise

        limits_s, limits = limits_fn(x)
        pha_df["limits_s"] = limits_s

        peaks_data: list[list[Any]] = []
        lineshapes = 'l'
        fited_curve: dict[str, np.ndarray] = {}
        error_curve: dict[str, np.ndarray] = {}
        rmse: dict[str, int] = {}

        # For each carbon/region, begin analysis
        for j, graph in enumerate(table):
            x_start = float(graph[0])
            x_end = float(graph[1])
            title = graph[2]

            row = j // 4
            col = j % 4

            # Find start and end in current data
            idx_start = None
            idx_end = None
            for index, el in enumerate(x):
                if min(x_start, x_end) <= el <= max(x_start, x_end):
                    if idx_start is None:
                        idx_start = index
                    idx_end = index

            x1: list[float] = []
            y1: list[float] = []
            if idx_start is not None and idx_end is not None:
                x1 = x[idx_start : idx_end + 1]
                y1 = y[idx_start : idx_end + 1]

            x1_np = np.array(x1)
            y1_np = np.array(y1)

            limits_s_peaks = [limit for limit in limits_s if limit[0] == title]
            for i, info in enumerate(limits_s_peaks):
                limits_s_peaks[i][2] = inv_scale(x1_np, x[int(info[2])])
                limits_s_peaks[i][3] = inv_scale(x1_np, x[int(info[3])])

            r = sum(1 for nb in y1_np if nb >= noise)

            peak_ppms: list[float] = []
            peak_intensities: list[float] = []

            if r >= 1:
                print(title, "Peaks have been found")
                peaks = ng.peakpick.pick(y1_np, pthres=noise, algorithm='downward', lineshapes=lineshapes)
                peak_index = peaks['X_AXIS'].astype(int)

                for index in peak_index:
                    peak_ppms.append(x1[index])
                    peak_intensities.append(y1[index])

                for bname, (end, start) in limits.items():
                    end_idx = inv_scale(x1_np, x[int(end)])
                    start_idx = inv_scale(x1_np, x[int(start)])
                    params = [[(px, max(0.1, lw))] for px, lw in zip(peaks['X_AXIS'], peaks['X_LW']) if start_idx < px < end_idx]
                    npeaks = len(params)
                    if npeaks == 0:
                        continue

                    bounds = [[[(None, None), (0, None)]] for _ in range(npeaks)]
                    amps = [y1_np[int(px)] for px in peaks['X_AXIS'] if start_idx < px < end_idx]
                    ampbounds = [(0, None) for _ in range(npeaks)]

                    _bb = np.zeros(len(y1_np))
                    _bb[start_idx:end_idx] = y1_np[start_idx:end_idx]

                    # Note: With error_flag=False, fit_NDregion returns a 3-tuple (param_best, amp_best, ier)
                    # where ier is an int status code from scipy.optimize.leastsq (1-4 = success, >4 = error).
                    # nmrglue lacks type stubs/overloads, causing Pylance to infer a 3-or-5 tuple union.
                    params_best, amp_best, iers = ng.linesh.fit_NDregion(  # type: ignore
                        _bb, lineshapes, params, amps, bounds, ampbounds, wmask=None, error_flag=False
                    )

                    if iers > 4:  # type: ignore
                        print(bname, "error", iers)

                    sdata = ng.linesh.sim_NDregion(_bb.shape, lineshapes, params_best, amp_best)
                    fited_curve[bname] = sdata
                    error_curve[bname] = sdata - _bb

                    region_err = error_curve[bname][start_idx:end_idx]
                    valid_err = region_err[~np.isnan(region_err)]
                    if len(valid_err) > 0:
                        mean_sq = np.square(valid_err).sum() / len(valid_err)
                        if not np.isnan(mean_sq) and mean_sq >= 0:
                            rmse[bname] = int(sqrt(mean_sq))
                        else:
                            rmse[bname] = 0
                    else:
                        rmse[bname] = 0

                    defined_regions = [limit for limit in limits_s if limit[0] == bname]
                    peak_y = []
                    for _, name, r_end, r_start in defined_regions:
                        for param, amp in zip(params_best, amp_best):
                            _peak = param[0][0]
                            fwhm = param[0][1]
                            if r_start < _peak < r_end and fwhm > 0:
                                peak_y.append([_peak, y1_np[int(_peak)], fwhm])
                            else:
                                peak_y.append([float(0), float(0), float(0)])

                        sorted_peak_y = sorted(peak_y, key=lambda px: px[1], reverse=True)
                        peaks_data.append([
                            bname,
                            name,
                            sorted_peak_y[0][0],
                            y1_np[int(sorted_peak_y[0][0])],
                            _bb[r_start:r_end].sum(),
                            pi / 2 * amp * sorted_peak_y[0][2],
                        ])
                        peak_y = []

            elif r == 0:
                for bname, (end, start) in limits.items():
                    if bname == title:
                        defined_regions = [limit for limit in limits_s if limit[0] == title]
                        for _, name, r_end, r_start in defined_regions:
                            peaks_data.append([bname, name, np.nan, np.nan, 0, 0])
                        peak_ppms = []
                        peak_intensities = []
                        error_curve[title] = -1 * y1_np
                        fited_curve[title] = np.zeros(len(y1_np))

            axes[row, col].plot(x1, y1, color="#004155", linewidth=1)

            if len(y1) > 0:
                y1_min = np.min(y1)
                y1_max = np.max(y1)
                margin_max = 0.5 * (y1_max - y1_min)
                margin_min = 0.1 * (y1_max - y1_min)
                axes[row, col].set_ylim(y1_min - margin_min, y1_max + margin_max)

            axes[row, col].set_xlim(min(x_start, x_end), max(x_start, x_end))
            axes[row, col].set_title(title)
            axes[row, col].grid(visible=True, color='0.5', alpha=0.5)
            axes[row, col].get_yaxis().set_visible(True)

            if title in error_curve:
                axes[row, col].plot(x1, error_curve[title], color="#ff6600", linewidth=0.5)
            if title in fited_curve:
                axes[row, col].fill(x1, fited_curve[title], color="#838383")
                axes[row, col].plot(x1, fited_curve[title], color="#91df00")
            axes[row, col].invert_xaxis()

            defined_regions = [limit for limit in limits_s if limit[0] == title]
            for _, name, r_start, r_end in defined_regions:
                r_start_ppm, r_end_ppm = x1[int(r_start)], x1[int(r_end)]
                height_peaks = []
                for pos_p, h_p in zip(peak_ppms, peak_intensities):
                    if min(r_start_ppm, r_end_ppm) < pos_p <= max(r_start_ppm, r_end_ppm):
                        height_peaks.append((h_p, pos_p))
                if height_peaks:
                    height_peaks_sorted = sorted(height_peaks, key=lambda px: px[0], reverse=True)
                    c = "#FF0000"
                    if monomer_colors:
                        for m_key, color in monomer_colors.items():
                            if f"{m_key}*" in name:
                                c = color
                                break
                    else:
                        if 'B*' in name:
                            c = "#0011ff"
                        else:
                            c = '#287743'

                    peak_height, peak_pos = height_peaks_sorted[0]
                    axes[row, col].annotate(
                        name,
                        xy=(peak_pos, peak_height),
                        xytext=(peak_pos, peak_height + 0.2 * peak_height),
                        ha='left',
                        va='top',
                        color=c,
                    )
                    axes[row, col].plot(peak_pos, peak_height, 'x', color=c)

        for i in range(len(table), 4 * n_rows):
            axes[i // 4, i % 4].axis('off')

        pha_df["peaks_data"] = peaks_data
        plt.tight_layout()

    except FileNotFoundError:
        print(f"Error: File {file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return x, y, pha_df, fig
