"""
Copolymer microstructure statistical modeling for PHACS-NMR.
Implements Kamiya (1989) Model 1 (Bernoullian), Model 2 (Markovian 1st order),
Model 3 (Copolymer blend via scipy.optimize.fsolve), Kamiya D parameter, and
Stanford ChemE 160 chi and number-average block length metrics.

References:
    Kamiya, N. et al. (1989). Macromolecules, 22(4), 1676-1682.
    Virtanen, P., Gommers, R. et al. (2020). SciPy 1.0: Fundamental Algorithms
    for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
"""
from typing import Any, cast

import numpy as np
import scipy.optimize as sp_opt


def calculate_d_parameter(diads: dict[str, float], aa_key: str, bb_key: str, ab_key: str, ba_key: str) -> float:
    """
    Kamiya 1989 Eq (2): D = (F_AA * F_BB) / (F_AB * F_BA).
    Returns 100.0 sentinel when division by zero occurs (no cross-diads).
    """
    ab = diads.get(ab_key, 0.0)
    ba = diads.get(ba_key, 0.0)
    if ab == 0 or ba == 0 or np.isnan(ab) or np.isnan(ba):
        return 100.0
    denom = ab * ba
    if denom == 0:
        return 100.0
    return float((diads.get(aa_key, 0.0) * diads.get(bb_key, 0.0)) / denom)


def calculate_chi_parameter(diads: dict[str, float], a_key: str, b_key: str, ba_key: str) -> float:
    """
    Stanford ChemE 160: chi = P2{BA} / (P1{A} * P1{B}).
    Returns NaN when monomer fraction is 0.
    """
    a = diads.get(a_key, 0.0)
    b = diads.get(b_key, 0.0)
    if a == 0 or b == 0 or np.isnan(a) or np.isnan(b):
        return float(np.nan)
    denom = a * b
    if denom == 0:
        return float(np.nan)
    return float(diads.get(ba_key, 0.0) / denom)


def calculate_average_lengths(
    diads: dict[str, float],
    a_key: str,
    b_key: str,
    ab_key: str,
    ba_key: str,
) -> tuple[float, float]:
    """
    Stanford ChemE 160:
    len_B = B* / B*A (or BA cross-diad)
    len_A = A* / A*B (or AB cross-diad)
    """
    ba = diads.get(ba_key, 0.0)
    ab = diads.get(ab_key, 0.0)
    b = diads.get(b_key, 0.0)
    a = diads.get(a_key, 0.0)

    len_b = float(b / ba) if (ba != 0 and not np.isnan(ba)) else float(np.nan)
    len_a = float(a / ab) if (ab != 0 and not np.isnan(ab)) else float(np.nan)
    return len_a, len_b


def compute_kamiya_model_1(avg_ratio: dict[str, float], keys_map: dict[str, str]) -> dict[str, float]:
    """
    Kamiya 1989 Model 1: Bernoullian random copolymer statistics.
    """
    f_a = avg_ratio[keys_map["A*"]]
    f_b = avg_ratio[keys_map["B*"]]

    m1 = {}
    m1[keys_map["AA"]] = f_a ** 2
    m1[keys_map["AB"]] = f_a * f_b
    m1[keys_map["BA"]] = m1[keys_map["AB"]]
    m1[keys_map["BB"]] = f_b ** 2

    m1[keys_map["AAA"]] = f_a ** 3
    m1[keys_map["BAA"]] = (f_a ** 2) * f_b
    m1[keys_map["AAB"]] = m1[keys_map["BAA"]]
    m1[keys_map["BAB"]] = f_a * (f_b ** 2)

    m1[keys_map["B*"]] = m1[keys_map["BB"]] + m1[keys_map["BA"]]
    m1[keys_map["A*"]] = m1[keys_map["AA"]] + m1[keys_map["AB"]]

    m1[keys_map["BBB"]] = f_b ** 3
    m1[keys_map["ABB"]] = f_a * (f_b ** 2)
    m1[keys_map["BBA"]] = m1[keys_map["ABB"]]
    m1[keys_map["ABA"]] = (f_a ** 2) * f_b

    return m1


def compute_kamiya_model_2(avg_ratio: dict[str, float], keys_map: dict[str, str]) -> dict[str, float]:
    """
    Kamiya 1989 Model 2: First-order Markovian copolymer statistics.
    """
    f_a = avg_ratio[keys_map["A*"]]
    f_b = avg_ratio[keys_map["B*"]]

    paa = avg_ratio[keys_map["AA"]] / f_a if f_a != 0 else 0.0
    pbb = avg_ratio[keys_map["BB"]] / f_b if f_b != 0 else 0.0
    pab = avg_ratio[keys_map["AB"]] / f_a if f_a != 0 else 0.0
    pba = avg_ratio[keys_map["BA"]] / f_b if f_b != 0 else 0.0

    denom = pab + pba
    if denom == 0:
        denom = 1e-12

    m2 = {}
    m2[keys_map["BB"]] = (pab * pbb) / denom
    m2[keys_map["BA"]] = (pab * pba) / denom
    m2[keys_map["AB"]] = m2[keys_map["BA"]]
    m2[keys_map["AA"]] = (paa * pba) / denom

    m2[keys_map["AAA"]] = (paa * paa * pba) / denom
    m2[keys_map["BAA"]] = (paa * pab * pba) / denom
    m2[keys_map["AAB"]] = m2[keys_map["BAA"]]
    m2[keys_map["BAB"]] = (pab * pab * pba) / denom

    m2[keys_map["B*"]] = m2[keys_map["BB"]] + m2[keys_map["BA"]]
    m2[keys_map["A*"]] = m2[keys_map["AA"]] + m2[keys_map["AB"]]

    m2[keys_map["BBB"]] = (pbb * pbb * pab) / denom
    m2[keys_map["ABB"]] = (pbb * pba * pab) / denom
    m2[keys_map["BBA"]] = m2[keys_map["ABB"]]
    m2[keys_map["ABA"]] = (pba * pba * pab) / denom

    m2[keys_map["Paa"]] = paa
    m2[keys_map["Pbb"]] = pbb
    m2[keys_map["Pab"]] = pab
    m2[keys_map["Pba"]] = pba

    return m2


def compute_kamiya_model_3(avg_ratio: dict[str, float], keys_map: dict[str, str]) -> dict[str, float]:
    """
    Kamiya 1989 Model 3: Copolymer blend solved via sp.optimize.fsolve.
    """
    target_aaa = avg_ratio[keys_map["AAA"]]
    target_baa = avg_ratio[keys_map["BAA"]]
    target_aab = avg_ratio[keys_map["AAB"]]
    target_bab = avg_ratio[keys_map["BAB"]]

    def functions_to_solve(z):
        a_val = z[0]
        b_val = z[1]
        x_val = z[2]
        f = np.empty(3)
        f[0] = (a_val ** 3) * x_val + (b_val ** 3) * (1 - x_val) - target_aaa
        f[1] = (a_val ** 2) * (1 - a_val) * x_val + (b_val ** 2) * (1 - b_val) * (1 - x_val) - (target_baa + target_aab) / 2
        f[2] = a_val * ((1 - a_val) ** 2) * x_val + b_val * ((1 - b_val) ** 2) * (1 - x_val) - target_bab
        return f

    z_guess = np.array([0, 0, 0])
    sol = cast(np.ndarray, sp_opt.fsolve(functions_to_solve, z_guess))
    a_sol, b_sol, x_sol = float(sol[0]), float(sol[1]), float(sol[2])


    m3 = {}
    m3[keys_map["AAA"]] = (a_sol ** 3) * x_sol + (b_sol ** 3) * (1 - x_sol)
    m3[keys_map["BAA"]] = (a_sol ** 2) * (1 - a_sol) * x_sol + (b_sol ** 2) * (1 - b_sol) * (1 - x_sol)
    m3[keys_map["AAB"]] = m3[keys_map["BAA"]]
    m3[keys_map["BAB"]] = a_sol * ((1 - a_sol) ** 2) * x_sol + b_sol * ((1 - b_sol) ** 2) * (1 - x_sol)

    m3[keys_map["BB"]] = ((1 - a_sol) ** 2) * x_sol + ((1 - b_sol) ** 2) * (1 - x_sol)
    m3[keys_map["BA"]] = a_sol * (1 - a_sol) * x_sol + b_sol * (1 - b_sol) * (1 - x_sol)
    m3[keys_map["AB"]] = m3[keys_map["BA"]]
    m3[keys_map["AA"]] = (a_sol ** 2) * x_sol + (b_sol ** 2) * (1 - x_sol)

    m3[keys_map["A*"]] = a_sol * x_sol + b_sol * (1 - x_sol)
    m3[keys_map["B*"]] = (1 - a_sol) * x_sol + (1 - b_sol) * (1 - x_sol)

    m3[keys_map["BBB"]] = ((1 - a_sol) ** 3) * x_sol + ((1 - b_sol) ** 3) * (1 - x_sol)
    m3[keys_map["ABB"]] = ((1 - a_sol) ** 2) * a_sol * x_sol + ((1 - b_sol) ** 2) * b_sol * (1 - x_sol)
    m3[keys_map["BBA"]] = m3[keys_map["ABB"]]
    m3[keys_map["ABA"]] = (1 - a_sol) * (a_sol ** 2) * x_sol + (1 - b_sol) * (b_sol ** 2) * (1 - x_sol)

    m3["A"] = a_sol
    m3["B"] = b_sol
    m3["X"] = x_sol

    return m3


def apply_copolymer_statistics(ratio_dic: dict[Any, dict[str, float]], keys_map: dict[str, str]) -> dict[Any, dict[str, float]]:
    """
    Compute Kamiya D, chi, len_A, len_B for each row (1, 2, 3, 4, avg, model_1, model_2, model_3).
    """
    aa_key = keys_map["AA"]
    bb_key = keys_map["BB"]
    ab_key = keys_map["AB"]
    ba_key = keys_map["BA"]
    a_key = keys_map["A*"]
    b_key = keys_map["B*"]
    len_a_key = keys_map["len_A"]
    len_b_key = keys_map["len_B"]

    for row_key, row in ratio_dic.items():
        row["D"] = calculate_d_parameter(row, aa_key, bb_key, ab_key, ba_key)
        row["chi"] = calculate_chi_parameter(row, a_key, b_key, ba_key)
        len_a, len_b = calculate_average_lengths(row, a_key, b_key, ab_key, ba_key)
        row[len_a_key] = len_a
        row[len_b_key] = len_b

    return ratio_dic
