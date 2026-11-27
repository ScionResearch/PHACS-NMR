"""
File I/O utilities for PHACS-NMR: Bruker convbin2asc ASCII spectrum files
and TopSpin deconvoluted (dconpeaks.txt) peak reports.
Formatted ASCII console tables rendered via tabulate (Sergey Astanin).
"""
import csv
import datetime
import os
from pathlib import Path
from typing import Any

from tabulate import tabulate

from phacs_nmr.path_utils import resolve_output_prefix


def read_ascii_spectrum(
    filename: str | Path,
    col_x: int = 3,
    col_y: int = 1,
    encoding: str = "cp1252",
) -> tuple[list[float], list[float]]:
    """
    Open the text file ascii-spec obtained by 'convbin2asc' command on TopSpin
    and convert it into two lists x (ppm) and y (intensity).
    """
    x: list[float] = []
    y: list[float] = []

    with open(str(filename), "r", encoding=encoding) as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split()
            if len(parts) == 4:
                cleaned = [part.replace(",", "") for part in parts]
                try:
                    x_val = float(cleaned[col_x])
                    y_val = float(cleaned[col_y])
                    x.append(x_val)
                    y.append(y_val)
                except ValueError:
                    continue

    return x, y


# Backwards compatibility alias
treatment_of_txt_file = read_ascii_spectrum


def read_dcon_file(filename: str | Path, encoding: str = "cp1252") -> dict[str, list[Any]]:
    """
    Read the text file containing the information of the deconvoluted spectrum (dconpeaks.txt).
    """
    with open(str(filename), "r", encoding=encoding) as f:
        lines = f.readlines()
        rows = [item.split() for item in lines[8:]]

        cleaned_rows = []
        for row in rows:
            if not row or (row and row[0] == "STD") or len(row) < 7:
                cleaned_rows.append([])
            else:
                cleaned_rows.append(row)

        peaks_ppm = []
        peaks_intensity = []
        peaks_area = []
        for row in cleaned_rows:
            if row:
                peaks_ppm.append(row[0])
                peaks_intensity.append(row[4])
                peaks_area.append(row[5])

        return {
            "peaks_ppm": peaks_ppm,
            "peaks_intensity": peaks_intensity,
            "peaks_area": peaks_area,
        }


# Backwards compatibility alias
data_dcon = read_dcon_file


def find_monomers(dconv_data: dict[str, list[Any]]) -> tuple[list[int], list[float], list[float]]:
    """
    Search for specific ppm values associated with known PHA monomers
    (F, B, V, X, Hp, O, N, >C9, >C9_1).
    """
    peaks_ppm = dconv_data["peaks_ppm"]
    peaks_intensity = dconv_data["peaks_intensity"]
    peaks_area = dconv_data["peaks_area"]

    limits = {
        "B": {"Est": (5.2608, 5.2631), "Met": (1.269, 1.2703), "CH2": (2.5965, 2.5985)},
        "V": {"Met": (0.8940, 0.897), "Est": (5.1516, 5.1554), "CH2": (2.556, 2.5595)},
        "X": {"Est": (5.204, 5.206), "Met": (0.9099, 0.9133), "CH2": (2.554, 2.556)},
        "Hp": {"Met": (0.9005, 0.9015)},
        "O": {"Est": (5.1809, 5.1834), "Met": (0.8778, 0.8802), "CH2": (2.550, 2.5474)},
        "N": {"Met": (0.8755, 0.8765)},
        "F": {"Est": (4.107, 4.110), "CH2": (2.350, 2.354)},
        ">C9": {"Met": (0.8764, 0.8774)},
        ">C9_1": {"Met": (0.881, 0.8825)},
    }

    keys = ["F", "B", "V", "X", "Hp", "O", "N", ">C9", ">C9_1"]
    p_id = [0] * 9
    est_monomers_int = [0.0] * 9
    est_monomers_area = [0.0] * 9

    for _id, position in enumerate(peaks_ppm):
        pos_f = float(position)
        for i, key in enumerate(keys):
            for proton, (mini, maxi) in limits[key].items():
                if proton == "Est" and min(float(mini), float(maxi)) <= pos_f <= max(float(mini), float(maxi)):
                        p_id[i] = 1
                        est_monomers_int[i] = float(peaks_intensity[_id])
                        est_monomers_area[i] = float(peaks_area[_id])

    total_int = sum(est_monomers_int)
    if total_int > 0:
        est_monomers_int[0] = float(est_monomers_int[0]) / total_int
    else:
        est_monomers_int[0] = 0.0

    print("estimated proportion of 4HB:", est_monomers_int[0])

    # First test: 4HB peak is found and intensity is high enough
    if p_id[0] == 1 and est_monomers_int[0] >= 0.01:
        print("The polymer is containing 4HB, the analysis is now done on the range [2ppm;3ppm]")
        for _id, position in enumerate(peaks_ppm):
            pos_f = float(position)
            for i, key in enumerate(keys):
                for proton, (mini, maxi) in limits[key].items():
                    if proton == "CH2" and min(float(mini), float(maxi)) <= pos_f <= max(float(mini), float(maxi)):
                            p_id[i] = 1
                            est_monomers_int[i] = float(peaks_intensity[_id])
                            est_monomers_area[i] = float(peaks_area[_id])

        if p_id == [1, 1, 0, 0, 0, 0, 0, 0, 0]:
            print("The polymer is containing only 3HB and 4HB, the analysis is now done on the range [4ppm;6ppm]")
            for _id, position in enumerate(peaks_ppm):
                pos_f = float(position)
                for i, key in enumerate(keys):
                    for proton, (mini, maxi) in limits[key].items():
                        if proton == "Est" and min(float(mini), float(maxi)) <= pos_f <= max(float(mini), float(maxi)):
                                p_id[i] = 1
                                est_monomers_int[i] = float(peaks_intensity[_id])
                                est_monomers_area[i] = float(peaks_area[_id])

    # Second test: 4HB peak is found but signal intensity is too low
    elif p_id[0] == 1 and est_monomers_int[0] < 0.01:
        p_id[0] = 0
        est_monomers_int[0] = 0.0
        print("4HB found but in very low proportion, work will be done on the range [0.5ppm-1.5ppm]")
        for _id, position in enumerate(peaks_ppm):
            pos_f = float(position)
            for i, key in enumerate(keys):
                for proton, (mini, maxi) in limits[key].items():
                    if proton == "Met" and min(float(mini), float(maxi)) <= pos_f <= max(float(mini), float(maxi)):
                            p_id[i] = 1
                            est_monomers_int[i] = float(peaks_intensity[_id])
                            est_monomers_area[i] = float(peaks_area[_id])

    # Third test: 4HB peak is not found
    elif p_id[0] == 0:
        print("4HB is not found, work will be done on the range [0.5ppm-1.5ppm]")
        for _id, position in enumerate(peaks_ppm):
            pos_f = float(position)
            for i, key in enumerate(keys):
                for proton, (mini, maxi) in limits[key].items():
                    if proton == "Met" and min(float(mini), float(maxi)) <= pos_f <= max(float(mini), float(maxi)):
                            p_id[i] = 1
                            est_monomers_int[i] = float(peaks_intensity[_id])
                            est_monomers_area[i] = float(peaks_area[_id])

    return p_id, est_monomers_area, est_monomers_int


def results_dconv(
    filename: str | Path,
    sensitivity: float = 0.01,
) -> tuple[list[str], list[int], str, list[float], list[float]]:
    """
    Estimate monomer composition from 1H NMR deconvoluted spectrum and export CSV report.
    """
    dconv_data = read_dcon_file(filename)
    p_id, est_monomers_area, est_monomers_int = find_monomers(dconv_data)
    monomers = ["F", "B", "V", "X", "Hp", "O", "N", ">C9", ">C9_1"]
    polymer_type = "PH"

    est_monomers_area = [float(v) for v in est_monomers_area]
    est_monomers_int = [float(v) for v in est_monomers_int]

    total_area = sum(est_monomers_area)
    total_intensities = sum(est_monomers_int)
    print("total_area", total_area)
    print("total_intensities", total_intensities)

    table_int = tabulate([est_monomers_int], headers=monomers, tablefmt="pretty", numalign="center")
    print("Intensities:", table_int)
    print(f"Sensibility of the analytical model: {sensitivity}")
    print("Calculation of the relative intensities...")

    if len(est_monomers_area) == len(est_monomers_int) == len(p_id):
        relative_total = []
        for i in range(len(p_id)):
            if p_id[i] == 1:
                est_monomers_int[i] = est_monomers_int[i] / total_intensities if total_intensities > 0 else 0.0
                est_monomers_area[i] = est_monomers_area[i] / total_area if total_area > 0 else 0.0

                if est_monomers_int[i] > float(sensitivity):
                    p_id[i] = 1
                    polymer_type = polymer_type + "-" + monomers[i]
                    relative_total.append(est_monomers_int[i])
                else:
                    p_id[i] = 0
                    est_monomers_int[i] = 0.0
                    relative_total.append(0.0)
            else:
                est_monomers_int[i] = 0.0
                est_monomers_area[i] = 0.0

        total_approx_int = sum(relative_total)
        for i in range(len(p_id)):
            est_monomers_int[i] = est_monomers_int[i] / total_approx_int if total_approx_int > 0 else 0.0

    headers = ["4HB", "3HB", "3HV", "3HHx", "3HHp", "3HO", "3HN", ">C9", ">C9_1"]

    if p_id != [0, 0, 0, 0, 0, 0, 0, 0, 0]:
        if p_id == [1, 1, 0, 0, 0, 0, 0, 0, 0]:
            est_monomers_int[0] = round(0.5437 * est_monomers_int[0], ndigits=3)
            est_monomers_int[1] = 1 - round(est_monomers_int[0], ndigits=3)
        elif p_id == [1, 1, 1, 0, 0, 0, 0, 0, 0]:
            est_monomers_int[0] = round(0.5437 * est_monomers_int[0], ndigits=3)
            est_monomers_int[2] = round(est_monomers_int[2], ndigits=3)
            est_monomers_int[1] = 1 - round(est_monomers_int[0] + est_monomers_int[2], ndigits=3)
        elif p_id == [0, 1, 1, 0, 0, 0, 0, 0, 0]:
            est_monomers_int[2] = round(0.8401 * est_monomers_int[2], ndigits=3)
            est_monomers_int[1] = 1 - round(est_monomers_int[2], ndigits=3)
        elif p_id == [0, 1, 0, 1, 0, 1, 0, 0, 0]:
            est_monomers_int[5] = round(1.1263 * est_monomers_int[5] + 0.0316, ndigits=3)
            est_monomers_int[1] = round(1.0344 * est_monomers_int[1] - 0.0686, ndigits=3)
            est_monomers_int[3] = round(1 - (est_monomers_int[1] + est_monomers_int[5]), ndigits=3)
        else:
            for i in range(len(est_monomers_int)):
                est_monomers_int[i] = round(est_monomers_int[i], ndigits=3)

            est_monomers_int[0] = round(0.5437 * est_monomers_int[0], ndigits=3)

            if est_monomers_int[1] != 0:
                est_monomers_int[1] = 1 - sum([
                    est_monomers_int[0], est_monomers_int[2], est_monomers_int[3],
                    est_monomers_int[4], est_monomers_int[5], est_monomers_int[6],
                    est_monomers_int[7], est_monomers_int[8],
                ])

        est_data = [est_monomers_int]

        print("Your polymer sample is a:", polymer_type)
        print(p_id, "\n", monomers)
        print("The estimation of the relative proportion of each monomer is obtained by analysing comparable peaks using dconv method in TopSpin")
        print("Calculation is using the the intensity of the mane peak for each hydrogen signal")
        print("Estimation results:")
        table = tabulate(
            est_data,
            headers=headers,
            tablefmt="pretty",
            numalign="center",
            colalign=("left", "center", "center", "center", "center", "center", "center"),
        )
        print(table)

        date = datetime.datetime.now()
        output_prefix = resolve_output_prefix(str(filename))
        output = f"{output_prefix}_H_s{sensitivity}_Software_results_{date.strftime('%d-%m-%Y %H-%M-%S')}.csv"
        final_path = os.path.join(os.path.dirname(str(filename)), output)

        with open(final_path, "w", newline="") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow([output])
            spamwriter.writerow("")
            spamwriter.writerow(["This excel sheet contains a estimation of the composition of your PHA sample."])
            spamwriter.writerow(["The estimation of the relative proportion of each monomer is obtained by analysing specific comparable peaks using dconv method in TopSpin"])
            spamwriter.writerow([f"The analysis is made at a level of {sensitivity} for the detection"])
            spamwriter.writerow(["F = 4-HydroxyButyrate"])
            spamwriter.writerow(["B = 3-HydroxyButyrate"])
            spamwriter.writerow(["V = 3-HydroxyValerate"])
            spamwriter.writerow(["X = 3-HydroxyHexanoate"])
            spamwriter.writerow(["Hp = 3-HydroxyHeptanoate"])
            spamwriter.writerow(["O = 3-HydroxyOctanoate"])
            spamwriter.writerow(["N = 3-HydroxyNonanoate"])
            spamwriter.writerow([">C9 = longer monomer than N"])
            spamwriter.writerow([">C9_1 = other longer monomer than N"])
            spamwriter.writerow([polymer_type])
            spamwriter.writerow(headers)
            spamwriter.writerow(est_monomers_int)
    else:
        print("No monomer has been found in the current list: \n", headers[1:])
        print("The polymer is an mcl-PHA composed of monomers with higher length than 3-HydroxyOctanoate")

    return headers, p_id, polymer_type, est_monomers_area, est_monomers_int
