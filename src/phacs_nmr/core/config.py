"""
Declarative configuration and PolymerProfile dataclasses for PHACS-NMR.
Defines spectral region tables, limit builders, resonance key mappings,
and plotting/export metadata for PHBV, P34HB, PHBHHx, and PHBHHxHO_MCL.
"""
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from phacs_nmr.core.noise import inv_scale

# =============================================================================
# Limit builders for each polymer type
# =============================================================================

def limits_phbv(ppm_scale: Any) -> tuple[list[list[Any]], dict[str, tuple[int, int]]]:
    """Define locations of individual resonances for PHBV."""
    limits_s = [
        ['V5', 'V*V', inv_scale(ppm_scale, 9.32), inv_scale(ppm_scale, 9.34)],
        ['V5', 'V*B', inv_scale(ppm_scale, 9.34), inv_scale(ppm_scale, 9.38)],
        ['B4', 'B*V', inv_scale(ppm_scale, 19.712), inv_scale(ppm_scale, 19.75)],
        ['B4', 'B*B', inv_scale(ppm_scale, 19.75), inv_scale(ppm_scale, 19.8)],
        ['V4', 'VV*V', inv_scale(ppm_scale, 26.74), inv_scale(ppm_scale, 26.775)],
        ['V4', 'BV*V', inv_scale(ppm_scale, 26.775), inv_scale(ppm_scale, 26.815)],
        ['V4', 'VV*B', inv_scale(ppm_scale, 26.815), inv_scale(ppm_scale, 26.846)],
        ['V4', 'BV*B', inv_scale(ppm_scale, 26.846), inv_scale(ppm_scale, 26.90)],
        ['V2', 'VV*V', inv_scale(ppm_scale, 38.639), inv_scale(ppm_scale, 38.67)],
        ['V2', 'BV*V', inv_scale(ppm_scale, 38.67), inv_scale(ppm_scale, 38.705)],
        ['V2', 'VV*B', inv_scale(ppm_scale, 38.764), inv_scale(ppm_scale, 38.795)],
        ['V2', 'BV*B', inv_scale(ppm_scale, 38.795), inv_scale(ppm_scale, 38.836)],
        ['B2', 'B*V', inv_scale(ppm_scale, 40.764), inv_scale(ppm_scale, 40.792)],
        ['B2', 'B*B', inv_scale(ppm_scale, 40.792), inv_scale(ppm_scale, 40.84)],
        ['B3', 'B*B', inv_scale(ppm_scale, 67.56), inv_scale(ppm_scale, 67.643)],
        ['B3', 'B*V', inv_scale(ppm_scale, 67.643), inv_scale(ppm_scale, 67.716)],
        ['V3', 'V*B', inv_scale(ppm_scale, 71.85), inv_scale(ppm_scale, 71.911)],
        ['V3', 'V*V', inv_scale(ppm_scale, 71.911), inv_scale(ppm_scale, 71.96)],
        ['B1V1', 'B*B', inv_scale(ppm_scale, 169.11), inv_scale(ppm_scale, 169.185)],
        ['B1V1', 'V*B', inv_scale(ppm_scale, 169.29), inv_scale(ppm_scale, 169.334)],
        ['B1V1', 'B*V', inv_scale(ppm_scale, 169.334), inv_scale(ppm_scale, 169.38)],
        ['B1V1', 'V*V', inv_scale(ppm_scale, 169.45), inv_scale(ppm_scale, 169.55)],
    ]

    limits = {
        'V2': (inv_scale(ppm_scale, 38.5), inv_scale(ppm_scale, 39)),
        'V4': (inv_scale(ppm_scale, 26.6), inv_scale(ppm_scale, 27)),
        'B4': (inv_scale(ppm_scale, 19.5), inv_scale(ppm_scale, 19.9)),
        'V5': (inv_scale(ppm_scale, 9.2), inv_scale(ppm_scale, 9.5)),
        'B1V1': (inv_scale(ppm_scale, 168.9), inv_scale(ppm_scale, 169.8)),
        'V3': (inv_scale(ppm_scale, 71.6), inv_scale(ppm_scale, 72.1)),
        'B3': (inv_scale(ppm_scale, 67.45), inv_scale(ppm_scale, 67.8)),
        'B2': (inv_scale(ppm_scale, 40.6), inv_scale(ppm_scale, 41)),
    }
    return limits_s, limits


def limits_p34hb(ppm_scale: Any) -> tuple[list[list[Any]], dict[str, tuple[int, int]]]:
    """Define locations of individual resonances for P34HB."""
    limits_s = [
        ['B4', 'BB*B', inv_scale(ppm_scale, 19.734), inv_scale(ppm_scale, 19.796)],
        ['B4', 'BB*F', inv_scale(ppm_scale, 19.796), inv_scale(ppm_scale, 19.83)],
        ['B4', 'FB*B*2', inv_scale(ppm_scale, 19.83), inv_scale(ppm_scale, 19.9)],
        ['B4', 'FB*F', inv_scale(ppm_scale, 19.9), inv_scale(ppm_scale, 19.977)],
        ['F3', 'BF*', inv_scale(ppm_scale, 23.9), inv_scale(ppm_scale, 23.985)],
        ['F3', 'FF*', inv_scale(ppm_scale, 23.985), inv_scale(ppm_scale, 24.05)],
        ['F2', 'BF*F', inv_scale(ppm_scale, 30.5), inv_scale(ppm_scale, 30.626)],
        ['F2', 'FF*F', inv_scale(ppm_scale, 30.626), inv_scale(ppm_scale, 30.684)],
        ['F2', 'BF*BB', inv_scale(ppm_scale, 30.75), inv_scale(ppm_scale, 30.81)],
        ['F2', 'BF*BF', inv_scale(ppm_scale, 30.81), inv_scale(ppm_scale, 30.841)],
        ['F2', 'FF*BB', inv_scale(ppm_scale, 30.841), inv_scale(ppm_scale, 30.868)],
        ['F2', 'FF*BF', inv_scale(ppm_scale, 30.868), inv_scale(ppm_scale, 30.91)],
        ['B2', 'BB*F', inv_scale(ppm_scale, 40.45), inv_scale(ppm_scale, 40.58)],
        ['B2', 'FB*F', inv_scale(ppm_scale, 40.61), inv_scale(ppm_scale, 40.73)],
        ['B2', 'BB*BB', inv_scale(ppm_scale, 40.76), inv_scale(ppm_scale, 40.84)],
        ['B2', 'BB*BF', inv_scale(ppm_scale, 40.84), inv_scale(ppm_scale, 40.885)],
        ['B2', 'FB*BB', inv_scale(ppm_scale, 40.885), inv_scale(ppm_scale, 40.935)],
        ['B2', 'FB*BF', inv_scale(ppm_scale, 40.935), inv_scale(ppm_scale, 41.1)],
        ['F4', 'FF*', inv_scale(ppm_scale, 63.47), inv_scale(ppm_scale, 63.56)],
        ['F4', 'BF*', inv_scale(ppm_scale, 63.56), inv_scale(ppm_scale, 63.68)],
        ['B3', 'FB*', inv_scale(ppm_scale, 67.38), inv_scale(ppm_scale, 67.46)],
        ['B3', 'BB*', inv_scale(ppm_scale, 67.56), inv_scale(ppm_scale, 67.7)],
        ['B1', 'BB*B', inv_scale(ppm_scale, 169.1), inv_scale(ppm_scale, 169.2)],
        ['B1', 'FB*B', inv_scale(ppm_scale, 169.2), inv_scale(ppm_scale, 169.3)],
        ['B1', 'BB*F', inv_scale(ppm_scale, 169.96), inv_scale(ppm_scale, 170.04)],
        ['B1', 'FB*F', inv_scale(ppm_scale, 170.04), inv_scale(ppm_scale, 170.12)],
        ['F1', 'F*B', inv_scale(ppm_scale, 171.8), inv_scale(ppm_scale, 171.93)],
        ['F1', 'BF*F', inv_scale(ppm_scale, 172.59), inv_scale(ppm_scale, 172.634)],
        ['F1', 'FF*F', inv_scale(ppm_scale, 172.634), inv_scale(ppm_scale, 172.7)],
    ]

    limits = {
        'B4': (inv_scale(ppm_scale, 19.65), inv_scale(ppm_scale, 20)),
        'F3': (inv_scale(ppm_scale, 23.8), inv_scale(ppm_scale, 24.2)),
        'F2': (inv_scale(ppm_scale, 30.5), inv_scale(ppm_scale, 31)),
        'B2': (inv_scale(ppm_scale, 40.45), inv_scale(ppm_scale, 41.2)),
        'F4': (inv_scale(ppm_scale, 63), inv_scale(ppm_scale, 64)),
        'B3': (inv_scale(ppm_scale, 67), inv_scale(ppm_scale, 68)),
        'B1': (inv_scale(ppm_scale, 168.5), inv_scale(ppm_scale, 170.5)),
        'F1': (inv_scale(ppm_scale, 171.5), inv_scale(ppm_scale, 173)),
    }
    return limits_s, limits


def limits_phbhhx(ppm_scale: Any) -> tuple[list[list[Any]], dict[str, tuple[int, int]]]:
    """Define locations of individual resonances for PHBHHx."""
    limits_s = [
        ['X2', 'X*X', inv_scale(ppm_scale, 39.1), inv_scale(ppm_scale, 39.25)],
        ['X2', 'X*B', inv_scale(ppm_scale, 39.25), inv_scale(ppm_scale, 39.5)],
        ['X4', 'BX*X', inv_scale(ppm_scale, 35.98), inv_scale(ppm_scale, 35.9995)],
        ['X4', 'XX*X', inv_scale(ppm_scale, 35.9995), inv_scale(ppm_scale, 36.01)],
        ['X4', 'XX*B', inv_scale(ppm_scale, 36.01), inv_scale(ppm_scale, 36.036)],
        ['X4', 'BX*B', inv_scale(ppm_scale, 36.036), inv_scale(ppm_scale, 36.054)],
        ['B4', 'XB*X', inv_scale(ppm_scale, 19.725), inv_scale(ppm_scale, 19.75)],
        ['B4', 'BB*X', inv_scale(ppm_scale, 19.75), inv_scale(ppm_scale, 19.7605)],
        ['B4', 'BB*B', inv_scale(ppm_scale, 19.7605), inv_scale(ppm_scale, 19.782)],
        ['B4', 'XB*B', inv_scale(ppm_scale, 19.782), inv_scale(ppm_scale, 19.82)],
        ['X5', 'X*X', inv_scale(ppm_scale, 18.32), inv_scale(ppm_scale, 18.336)],
        ['X5', 'X*B', inv_scale(ppm_scale, 18.336), inv_scale(ppm_scale, 18.37)],
        ['X6', 'X*', inv_scale(ppm_scale, 13.78), inv_scale(ppm_scale, 13.83)],
        ['B1X1', 'X*X', inv_scale(ppm_scale, 169.46), inv_scale(ppm_scale, 169.5)],
        ['B1X1', 'B*X', inv_scale(ppm_scale, 169.31), inv_scale(ppm_scale, 169.4)],
        ['B1X1', 'X*B', inv_scale(ppm_scale, 169.26), inv_scale(ppm_scale, 169.31)],
        ['B1X1', 'B*B', inv_scale(ppm_scale, 169.1), inv_scale(ppm_scale, 169.18)],
        ['X3', 'X*X', inv_scale(ppm_scale, 70.4), inv_scale(ppm_scale, 70.55)],
        ['X3', 'X*B', inv_scale(ppm_scale, 70.55), inv_scale(ppm_scale, 70.62)],
        ['B3', 'XB*B', inv_scale(ppm_scale, 67.56), inv_scale(ppm_scale, 67.61)],
        ['B3', 'BB*B', inv_scale(ppm_scale, 67.61), inv_scale(ppm_scale, 67.655)],
        ['B3', 'BB*X', inv_scale(ppm_scale, 67.655), inv_scale(ppm_scale, 67.678)],
        ['B3', 'XB*X', inv_scale(ppm_scale, 67.678), inv_scale(ppm_scale, 67.71)],
        ['B2', 'B*X', inv_scale(ppm_scale, 40.73), inv_scale(ppm_scale, 40.774)],
        ['B2', 'B*B', inv_scale(ppm_scale, 40.774), inv_scale(ppm_scale, 40.84)],
    ]

    limits = {
        'X2': (inv_scale(ppm_scale, 38.5), inv_scale(ppm_scale, 40)),
        'X4': (inv_scale(ppm_scale, 35), inv_scale(ppm_scale, 36.5)),
        'B4': (inv_scale(ppm_scale, 19.5), inv_scale(ppm_scale, 19.9)),
        'X5': (inv_scale(ppm_scale, 18.3), inv_scale(ppm_scale, 18.4)),
        'X6': (inv_scale(ppm_scale, 13.75), inv_scale(ppm_scale, 13.85)),
        'B1X1': (inv_scale(ppm_scale, 169), inv_scale(ppm_scale, 169.8)),
        'X3': (inv_scale(ppm_scale, 70.35), inv_scale(ppm_scale, 70.65)),
        'B3': (inv_scale(ppm_scale, 67.45), inv_scale(ppm_scale, 67.8)),
        'B2': (inv_scale(ppm_scale, 40.6), inv_scale(ppm_scale, 41)),
    }
    return limits_s, limits


def limits_phbhhxho_mcl(ppm_scale: Any) -> tuple[list[list[Any]], dict[str, tuple[int, int]]]:
    """Define locations of individual resonances for PHBHHxHO_MCL."""
    limits_s = [
        ['X6', 'X*', inv_scale(ppm_scale, 13.78), inv_scale(ppm_scale, 13.83)],
        ['O8', 'O*', inv_scale(ppm_scale, 13.935), inv_scale(ppm_scale, 14.093)],
        ['mcl', 'mcl*', inv_scale(ppm_scale, 14.093), inv_scale(ppm_scale, 14.16)],
        ['X5', 'X*X', inv_scale(ppm_scale, 18.32), inv_scale(ppm_scale, 18.336)],
        ['X5', 'X*B', inv_scale(ppm_scale, 18.336), inv_scale(ppm_scale, 18.37)],
        ['X5', 'X*O', inv_scale(ppm_scale, 18.41), inv_scale(ppm_scale, 18.47)],
        ['B4', 'B*', inv_scale(ppm_scale, 19.73), inv_scale(ppm_scale, 19.82)],
        ['O7', 'O*', inv_scale(ppm_scale, 22.47), inv_scale(ppm_scale, 22.66)],
        ['O6', 'O*O', inv_scale(ppm_scale, 24.66), inv_scale(ppm_scale, 24.745)],
        ['O6', 'O*B', inv_scale(ppm_scale, 24.745), inv_scale(ppm_scale, 24.8)],
        ['O5', 'O*O', inv_scale(ppm_scale, 31.48), inv_scale(ppm_scale, 31.62)],
        ['O5', 'O*B', inv_scale(ppm_scale, 31.62), inv_scale(ppm_scale, 31.68)],
        ['O4', 'O*B', inv_scale(ppm_scale, 33.69), inv_scale(ppm_scale, 33.73)],
        ['O4', 'O*O', inv_scale(ppm_scale, 33.73), inv_scale(ppm_scale, 33.777)],
        ['O4', 'O*X', inv_scale(ppm_scale, 33.777), inv_scale(ppm_scale, 33.81)],
        ['X4', 'X*O', inv_scale(ppm_scale, 35.86), inv_scale(ppm_scale, 35.95)],
        ['X4', 'X*X', inv_scale(ppm_scale, 35.98), inv_scale(ppm_scale, 36.01)],
        ['X4', 'X*B', inv_scale(ppm_scale, 36.01), inv_scale(ppm_scale, 36.054)],
        ['O2X2', 'O*', inv_scale(ppm_scale, 39), inv_scale(ppm_scale, 39.15)],
        ['O2X2', 'X*X', inv_scale(ppm_scale, 39.15), inv_scale(ppm_scale, 39.25)],
        ['O2X2', 'X*B', inv_scale(ppm_scale, 39.25), inv_scale(ppm_scale, 39.5)],
        ['B2', 'B*X', inv_scale(ppm_scale, 40.73), inv_scale(ppm_scale, 40.77)],
        ['B2', 'B*B', inv_scale(ppm_scale, 40.77), inv_scale(ppm_scale, 41)],
        ['B3', 'B*B', inv_scale(ppm_scale, 67.56), inv_scale(ppm_scale, 67.655)],
        ['B3', 'B*X', inv_scale(ppm_scale, 67.655), inv_scale(ppm_scale, 67.71)],
        ['X3', 'X*X', inv_scale(ppm_scale, 70.4), inv_scale(ppm_scale, 70.55)],
        ['X3', 'X*B', inv_scale(ppm_scale, 70.55), inv_scale(ppm_scale, 70.602)],
        ['X3', 'X*O', inv_scale(ppm_scale, 70.602), inv_scale(ppm_scale, 70.66)],
        ['O3', 'O*O*1', inv_scale(ppm_scale, 70.75), inv_scale(ppm_scale, 70.84)],
        ['O3', 'O*O*2', inv_scale(ppm_scale, 70.84), inv_scale(ppm_scale, 70.93)],
        ['B1X1O1', 'O*O', inv_scale(ppm_scale, 169.355), inv_scale(ppm_scale, 169.435)],
        ['B1X1O1', 'X*X', inv_scale(ppm_scale, 169.475), inv_scale(ppm_scale, 169.515)],
        ['B1X1O1', 'B*X', inv_scale(ppm_scale, 169.31), inv_scale(ppm_scale, 169.355)],
        ['B1X1O1', 'X*B', inv_scale(ppm_scale, 169.2), inv_scale(ppm_scale, 169.31)],
        ['B1X1O1', 'B*B', inv_scale(ppm_scale, 169), inv_scale(ppm_scale, 169.2)],
    ]

    limits = {
        'X6': (inv_scale(ppm_scale, 13.7), inv_scale(ppm_scale, 13.9)),
        'O8': (inv_scale(ppm_scale, 13.9), inv_scale(ppm_scale, 14.1)),
        'mcl': (inv_scale(ppm_scale, 14.093), inv_scale(ppm_scale, 14.16)),
        'X5': (inv_scale(ppm_scale, 18.2), inv_scale(ppm_scale, 18.6)),
        'B4': (inv_scale(ppm_scale, 19.6), inv_scale(ppm_scale, 19.9)),
        'O7': (inv_scale(ppm_scale, 22.4), inv_scale(ppm_scale, 22.7)),
        'O6': (inv_scale(ppm_scale, 24.6), inv_scale(ppm_scale, 24.9)),
        'O5': (inv_scale(ppm_scale, 30.6), inv_scale(ppm_scale, 31.8)),
        'O4': (inv_scale(ppm_scale, 33.65), inv_scale(ppm_scale, 33.85)),
        'X4': (inv_scale(ppm_scale, 34), inv_scale(ppm_scale, 37)),
        'O2X2': (inv_scale(ppm_scale, 38.5), inv_scale(ppm_scale, 39.5)),
        'B2': (inv_scale(ppm_scale, 40), inv_scale(ppm_scale, 42)),
        'B3': (inv_scale(ppm_scale, 67.4), inv_scale(ppm_scale, 67.8)),
        'X3': (inv_scale(ppm_scale, 70.5), inv_scale(ppm_scale, 70.7)),
        'O3': (inv_scale(ppm_scale, 70.7), inv_scale(ppm_scale, 71.0)),
        'B1X1O1': (inv_scale(ppm_scale, 169), inv_scale(ppm_scale, 169.6)),
    }
    return limits_s, limits


# =============================================================================
# PolymerProfile Dataclass
# =============================================================================

@dataclass
class PolymerProfile:
    """Declarative specification for a PHA copolymer analysis profile."""
    name: str
    display_name: str
    spectral_regions: list[tuple[float, float, str]]
    limits_builder: Callable[[Any], tuple[list[list[Any]], dict[str, tuple[int, int]]]]
    monomers: list[str]
    monomer_names: dict[str, str]
    monomer_colors: dict[str, str]
    keys_map: dict[str, str]
    csv_header: list[str] = field(default_factory=list)
    csv_notes: list[str] = field(default_factory=list)
    figsize_regions: tuple[float, float] = (10.0, 7.0)
    noise_as_int: bool = False
    has_mcl: bool = False
    has_csv: bool = True
    active_row_ids: list[Any] = field(default_factory=lambda: [1, 2, 3, 4, 'avg'])


# =============================================================================
# Standard Profiles Registration
# =============================================================================

PHBV_PROFILE = PolymerProfile(
    name="PHBV",
    display_name="Poly(3-hydroxybutyrate-co-3-hydroxyvalerate)",
    spectral_regions=[
        (168.9, 169.8, 'B1V1'),
        (71.6, 72.1, 'V3'),
        (67.0, 68.0, 'B3'),
        (40.3, 41.2, 'B2'),
        (38.5, 39.0, 'V2'),
        (26.6, 27.0, 'V4'),
        (19.65, 20.0, 'B4'),
        (9.2, 9.5, 'V5'),
    ],
    limits_builder=limits_phbv,
    monomers=['B', 'V'],
    monomer_names={'B': '3-HydroxyButyrate', 'V': '3-HydroxyValerate'},
    monomer_colors={'B': '#0011ff', 'V': '#287743'},
    keys_map={
        'A*': 'V*', 'B*': 'B*',
        'AA': 'V*V', 'AB': 'V*B', 'BA': 'B*V', 'BB': 'B*B',
        'AAA': 'VV*V', 'BAA': 'BV*V', 'AAB': 'VV*B', 'BAB': 'BV*B',
        'BBB': 'BB*B', 'ABB': 'VB*B', 'BBA': 'BB*V', 'ABA': 'VB*V',
        'Paa': 'Pvv', 'Pbb': 'Pbb', 'Pab': 'Pvb', 'Pba': 'Pbv',
        'len_A': 'len_V', 'len_B': 'len_B',
    },
    csv_header=[
        'id', 'D', 'chi', 'len_B', 'len_V',
        '', 'B*', 'B*B', 'B*V', 'V*', 'V*B', 'V*V',
        '', 'VV*V', 'BV*V', 'VV*B', 'BV*B', 'VB*V', 'BB*V', 'VB*B', 'BB*B',
        '', 'Pvv', 'Pbb', 'Pvb', 'Pbv',
        '', 'A', 'B', 'X'
    ],
    csv_notes=[
        'D = 1: Completely random copolymer',
        'D > 1: Block copolymer or mix of different random copolymer',
        'D < 1: Copolymer with an alternative tendency',
        'chi = 1: Completely random copolymer',
        'chi > 1: Copolymer with an alternating tendency',
        'chi = 2: Completely alternating copolymer',
        'chi < 1: Copolymer with a blocking tendency',
        'chi = 0: Completely block copolymer',
        'nan: could not calculate - division by 0',
        '',
        'Sources:',
        'Models and D values from: Microstructure of bacterially synthesized '
        'poly(3-hydroxybutyrate-co-3-hydroxyvalerate) - Naoko Kamiya - 1989',
        'Chi and average length of A and B runs as described in '
        'http://web.stanford.edu/class/cheme160/lectures/lecture7.pdf',
    ],
    figsize_regions=(8.27, 6.0),
    noise_as_int=True,
    has_mcl=False,
    has_csv=True,
    active_row_ids=[1, 2, 3, 4, 'avg'],
)

P34HB_PROFILE = PolymerProfile(
    name="P34HB",
    display_name="Poly(3-hydroxybutyrate-co-4-hydroxybutyrate)",
    spectral_regions=[
        (171.5, 173.0, 'F1'),
        (168.5, 170.5, 'B1'),
        (67.0, 68.0, 'B3'),
        (63.0, 64.0, 'F4'),
        (40.3, 41.2, 'B2'),
        (30.5, 31.0, 'F2'),
        (23.8, 24.2, 'F3'),
        (19.65, 20.0, 'B4'),
    ],
    limits_builder=limits_p34hb,
    monomers=['B', 'F'],
    monomer_names={'B': '3-HydroxyButyrate', 'F': '4-HydroxyButyrate'},
    monomer_colors={'B': '#0011ff', 'F': '#287743'},
    keys_map={
        'A*': 'F*', 'B*': 'B*',
        'AA': 'FF*', 'AB': 'FB*', 'BA': 'BF*', 'BB': 'BB*',
        'AAA': 'FF*F', 'BAA': 'BF*F', 'AAB': 'FF*B', 'BAB': 'BF*B',
        'BBB': 'BB*B', 'ABB': 'FB*B', 'BBA': 'BB*F', 'ABA': 'FB*F',
        'Paa': 'Pff', 'Pbb': 'Pbb', 'Pab': 'Pfb', 'Pba': 'Pbf',
        'len_A': 'len_F', 'len_B': 'len_B',
    },
    csv_header=[
        'id', 'D', 'chi', 'len_B', 'len_F',
        '', 'B*', 'BB*', 'FB*', 'F*', 'BF*', 'FF*',
        '', 'FF*F', 'BF*F', 'FF*B', 'BF*B', 'FB*F', 'BB*F', 'FB*B', 'BB*B',
        '', 'Pff', 'Pbb', 'Pfb', 'Pbf',
        '', 'A', 'B', 'X'
    ],
    csv_notes=[
        'chi = 1: Completely random copolymer',
        'chi > 1: Copolymer with an alternating tendency',
        'chi = 2: Completely alternating copolymer',
        'chi < 1: Copolymer with a blocking tendency',
        'chi = 0: Completely block copolymer',
        'nan: could not calculate - division by 0',
        '',
        'Sources:',
        'Models and D values from: Microstructure of bacterially synthesized '
        'poly(3-hydroxybutyrate-co-3-hydroxyvalerate) - Naoko Kamiya - 1989',
        'Chi and average length of A and B runs as described in '
        'http://web.stanford.edu/class/cheme160/lectures/lecture7.pdf',
    ],
    figsize_regions=(10.0, 7.0),
    noise_as_int=False,
    has_mcl=False,
    has_csv=True,
    active_row_ids=[1, 2, 3, 4, 'avg'],
)

PHBHHX_PROFILE = PolymerProfile(
    name="PHBHHx",
    display_name="Poly(3-hydroxybutyrate-co-3-hydroxyhexanoate)",
    spectral_regions=[
        (169.0, 169.8, 'B1X1'),
        (70.5, 71.0, 'X3'),
        (67.45, 67.8, 'B3'),
        (40.6, 41.0, 'B2'),
        (38.5, 40.0, 'X2'),
        (35.0, 36.5, 'X4'),
        (19.5, 19.9, 'B4'),
        (18.3, 18.4, 'X5'),
        (13.75, 13.85, 'X6'),
    ],
    limits_builder=limits_phbhhx,
    monomers=['B', 'X'],
    monomer_names={'B': '3-HydroxyButyrate', 'X': '3-HydroxyHexanoate'},
    monomer_colors={'B': '#0011ff', 'X': '#287743'},
    keys_map={
        'A*': 'X*', 'B*': 'B*',
        'AA': 'X*X', 'AB': 'X*B', 'BA': 'B*X', 'BB': 'B*B',
        'AAA': 'XX*X', 'BAA': 'BX*X', 'AAB': 'XX*B', 'BAB': 'BX*B',
        'BBB': 'BB*B', 'ABB': 'XB*B', 'BBA': 'BB*X', 'ABA': 'XB*X',
        'Paa': 'Pxx', 'Pbb': 'Pbb', 'Pab': 'Pxb', 'Pba': 'Pbx',
        'len_A': 'len_X', 'len_B': 'len_B',
    },
    csv_header=[
        'id', 'D', 'chi', 'len_B', 'len_X',
        '', 'B*', 'B*B', 'B*X', 'X*', 'X*B', 'X*X',
        '', 'XX*X', 'BX*X', 'XX*B', 'BX*B', 'XB*X', 'BB*X', 'XB*B', 'BB*B',
        '', 'Pxx', 'Pbb', 'Pxb', 'Pbx',
        '', 'A', 'B', 'X'
    ],
    csv_notes=[
        'chi = 1: Completely random copolymer',
        'chi > 1: Copolymer with an alternating tendency',
        'chi = 2: Completely alternating copolymer',
        'chi < 1: Copolymer with a blocking tendency',
        'chi = 0: Completely block copolymer',
        'nan: could not calculate - division by 0',
        '',
        'Sources:',
        'Models and D values from: Microstructure of bacterially synthesized '
        'poly(3-hydroxybutyrate-co-3-hydroxyvalerate) - Naoko Kamiya - 1989',
        'Chi and average length of A and B runs as described in '
        'http://web.stanford.edu/class/cheme160/lectures/lecture7.pdf',
    ],
    figsize_regions=(10.0, 7.0),
    noise_as_int=False,
    has_mcl=False,
    has_csv=True,
    active_row_ids=[1, 2, 3, 4, 'avg'],
)

PHBHHXHO_MCL_PROFILE = PolymerProfile(
    name="PHBHHxHO_MCL",
    display_name="Poly(3-hydroxybutyrate-co-3-hydroxyhexanoate-co-3-hydroxyoctanoate-co-mcl)",
    spectral_regions=[
        (169.6, 169.0, 'B1X1O1'),
        (71.0, 70.7, 'O3'),
        (70.7, 70.5, 'X3'),
        (67.8, 67.4, 'B3'),
        (42.0, 40.0, 'B2'),
        (39.5, 38.5, 'O2X2'),
        (36.5, 35.5, 'X4'),
        (33.85, 33.65, 'O4'),
        (31.8, 30.6, 'O5'),
        (24.9, 24.6, 'O6'),
        (22.7, 22.4, 'O7'),
        (19.9, 19.6, 'B4'),
        (18.6, 18.2, 'X5'),
        (14.16, 14.093, 'mcl'),
        (14.1, 13.9, 'O8'),
        (13.9, 13.7, 'X6'),
    ],
    limits_builder=limits_phbhhxho_mcl,
    monomers=['B', 'X', 'O', 'mcl'],
    monomer_names={
        'B': '3-HydroxyButyrate',
        'X': '3-HydroxyHexanoate',
        'O': '3-HydroxyOctanoate',
        'mcl': 'Medium Chain Length (mcl)',
    },
    monomer_colors={
        'B': '#0011ff',
        'X': '#287743',
        'O': '#b35900',
        'mcl': '#990099',
    },
    keys_map={},
    figsize_regions=(10.0, 7.0),
    noise_as_int=False,
    has_mcl=True,
    has_csv=False,
)

PROFILES: dict[str, PolymerProfile] = {
    "PHBV": PHBV_PROFILE,
    "P34HB": P34HB_PROFILE,
    "PHBHHx": PHBHHX_PROFILE,
    "PHBHHxHO_mcl": PHBHHXHO_MCL_PROFILE,
    "PHBHHxHO_MCL": PHBHHXHO_MCL_PROFILE,
}


def get_profile(name_or_id: str) -> PolymerProfile:
    """Retrieve profile by polymer name."""
    norm = name_or_id.strip()
    if norm in PROFILES:
        return PROFILES[norm]
    for key, p in PROFILES.items():
        if key.lower() == norm.lower():
            return p
    raise KeyError(f"Unknown polymer profile: {name_or_id}")
