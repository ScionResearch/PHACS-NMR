"""
Unified analysis engine for PHACS-NMR.
Orchestrates spectral parsing, peak picking, lineshape fitting,
copolymer sequence modeling (Kamiya models 1, 2, 3), figure generation,
PDF reporting, and CSV export.

Uses NumPy, SciPy, pandas, Matplotlib, seaborn, ReportLab, and nmrglue.
Please cite Virtanen et al. (2020), Harris et al. (2020), Hunter (2007),
and Helmus & Jaroniec (2013) when publishing results from this engine.
"""
import csv
import datetime
import io
import os
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.platypus import Image

from phacs_nmr.core.config import PolymerProfile
from phacs_nmr.core.models import (
    apply_copolymer_statistics,
    compute_kamiya_model_1,
    compute_kamiya_model_2,
    compute_kamiya_model_3,
)
from phacs_nmr.core.peak_fitting import split_graphic_regions
from phacs_nmr.path_utils import resolve_output_prefix, resolve_sample_name

# =============================================================================
# Copolymer Ratio & Model Computations
# =============================================================================

def compute_area_and_ratio(pha_df: dict[str, Any], profile: PolymerProfile, use_line_fitting: bool = False) -> dict[str, Any]:
    """Compute monomer and sequence areas and ratios based on raw or fitted peaks."""
    peaks_data = pha_df.get('peaks_data', [])
    area_used = 'Area_curve' if use_line_fitting else 'Area_raw'

    if not peaks_data:
        pha_df['peaks_data_df'] = pd.DataFrame(
            columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve']
        )
        pha_df['Ratio_dic'] = {}
        if profile.name == "PHBV":
            pha_df['Ratio_dic']['avg'] = {
                'D': 0, 'B*': 0, 'V*': 0, 'B*B': 0, 'B*V': 0, 'V*B': 0, 'V*V': 0,
                'VV*V': 0, 'BV*V': 0, 'VV*B': 0, 'BV*B': 0, 'chi': 0, 'len_B': 0, 'len_V': 0
            }
        elif profile.name == "P34HB":
            pha_df['Ratio_dic']['avg'] = {
                'D': 0, 'B*': 0, 'F*': 0, 'BB*': 0, 'BF*': 0, 'FB*': 0, 'FF*': 0,
                'FF*F': 0, 'BF*F': 0, 'FF*B': 0, 'BF*B': 0, 'BB*F': 0, 'FB*F': 0,
                'BB*B': 0, 'FB*B': 0, 'chi': 0, 'len_B': 0, 'len_F': 0
            }
        elif profile.name == "PHBHHx":
            pha_df['Ratio_dic']['avg'] = {
                'D': 0, 'B*': 0, 'X*': 0, 'B*B': 0, 'B*X': 0, 'X*B': 0, 'X*X': 0,
                'XX*X': 0, 'BX*X': 0, 'XX*B': 0, 'BX*B': 0, 'BB*B': 0, 'XB*B': 0,
                'BB*X': 0, 'XB*X': 0, 'chi': 0, 'len_B': 0, 'len_X': 0
            }
        return pha_df

    df = pd.DataFrame(peaks_data, columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve'])
    Ratio_dic: dict[Any, dict[str, float]] = {}

    if profile.name == "PHBV":
        df['MainGroup'] = pd.Categorical(df['MainGroup'], ['B1V1', 'V3', 'B3', 'B2', 'V2', 'V4', 'B4', 'V5'])
        sf = df[['MainGroup', 'Name', area_used]]
        Area_dic: dict[str, pd.DataFrame] = {}
        Area_dic['b1v1'] = sf[sf.MainGroup == 'B1V1'].copy()
        sum_b1v1 = Area_dic['b1v1'][area_used].sum()
        Area_dic['b1v1']['perc'] = Area_dic['b1v1'][area_used] / sum_b1v1 if sum_b1v1 != 0 else 0.0

        keys = ['V*', 'B*', 'V*V', 'V*B', 'B*V', 'B*B', 'VV*V', 'BV*V', 'VV*B', 'BV*B']
        str_keys = ['V\*', 'B\*', 'V\*V', 'V\*B', 'B\*V', 'B\*B', 'VV\*V', 'BV\*V', 'VV\*B', 'BV\*B']

        ratio: dict[str, float] = {}
        for _idx, _key in enumerate(keys):
            _value = np.nan
            if _key in ('V*V', 'B*V', 'V*B', 'B*B', 'B*', 'V*'):
                _value = Area_dic['b1v1'][Area_dic['b1v1']['Name'].str.contains(str_keys[_idx])]['perc'].sum()
            ratio[_key] = _value
        Ratio_dic[1] = ratio

        for _id in range(2, 5):
            _b = 'B' + str(_id)
            _v = 'V' + str(_id)
            Area_dic[_b] = sf[sf.MainGroup == _b].copy()
            Area_dic[_v] = sf[sf.MainGroup == _v].copy()
            _total = Area_dic[_b][area_used].sum() + Area_dic[_v][area_used].sum()
            Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
            Area_dic[_v]['perc'] = Area_dic[_v][area_used] / _total if _total != 0 else 0.0

            _ratio = {
                'V*': Area_dic[_v][Area_dic[_v]['Name'].str.contains('V\*')]['perc'].sum(),
                'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
                'V*V': Area_dic[_v][Area_dic[_v]['Name'].str.contains('V\*V')]['perc'].sum(),
                'V*B': Area_dic[_v][Area_dic[_v]['Name'].str.contains('V\*B')]['perc'].sum(),
                'B*V': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*V')]['perc'].sum(),
                'B*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*B')]['perc'].sum(),
                'VV*V': np.nan, 'BV*V': np.nan, 'VV*B': np.nan, 'BV*B': np.nan,
            }
            if _id in (2, 4):
                _ratio['VV*V'] = Area_dic[_v][Area_dic[_v]['Name'].str.contains('VV\*V')]['perc'].sum()
                _ratio['BV*V'] = Area_dic[_v][Area_dic[_v]['Name'].str.contains('BV\*V')]['perc'].sum()
                _ratio['VV*B'] = Area_dic[_v][Area_dic[_v]['Name'].str.contains('VV\*B')]['perc'].sum()
                _ratio['BV*B'] = Area_dic[_v][Area_dic[_v]['Name'].str.contains('BV\*B')]['perc'].sum()
            Ratio_dic[_id] = _ratio

        Ratio_dic['avg'] = {}
        for _key in keys:
            Ratio_dic['avg'][_key] = float(np.mean([
                Ratio_dic[1][_key], Ratio_dic[2][_key], Ratio_dic[3][_key], Ratio_dic[4][_key]
            ]))
        for _key in ['VV*V', 'BV*V', 'VV*B', 'BV*B']:
            Ratio_dic['avg'][_key] = float(np.mean([Ratio_dic[2][_key], Ratio_dic[4][_key]]))

    elif profile.name == "P34HB":
        df['MainGroup'] = pd.Categorical(df['MainGroup'], ['B1', 'F1', 'B2', 'F2', 'B3', 'F4', 'B4', 'F3'])
        sf = df[['MainGroup', 'Name', area_used]]
        Area_dic: dict[str, pd.DataFrame] = {}
        keys = ['F*', 'B*', 'FF*', 'FB*', 'BF*', 'BB*', 'FF*F', 'BF*F', 'FF*B', 'BF*B', 'BB*F', 'FB*F', 'BB*B', 'FB*B']

        for _id in range(1, 3):
            _b = 'B' + str(_id)
            _f = 'F' + str(_id)
            Area_dic[_b] = sf[sf.MainGroup == _b].copy()
            Area_dic[_f] = sf[sf.MainGroup == _f].copy()
            _total = Area_dic[_b][area_used].sum() + Area_dic[_f][area_used].sum()
            Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
            Area_dic[_f]['perc'] = Area_dic[_f][area_used] / _total if _total != 0 else 0.0

            _ratio = {
                'F*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('F\*')]['perc'].sum(),
                'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
                'FF*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('F\*F')]['perc'].sum(),
                'FB*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('F\*B')]['perc'].sum(),
                'BF*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*F')]['perc'].sum(),
                'BB*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*B')]['perc'].sum(),
                'FF*F': Area_dic[_f][Area_dic[_f]['Name'].str.contains('FF\*F')]['perc'].sum(),
                'BF*F': Area_dic[_f][Area_dic[_f]['Name'].str.contains('BF\*F')]['perc'].sum(),
                'FF*B': Area_dic[_f][Area_dic[_f]['Name'].str.contains('FF\*B')]['perc'].sum(),
                'BF*B': Area_dic[_f][Area_dic[_f]['Name'].str.contains('BF\*B')]['perc'].sum(),
                'BB*F': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*F')]['perc'].sum(),
                'FB*F': Area_dic[_b][Area_dic[_b]['Name'].str.contains('FB\*F')]['perc'].sum(),
                'BB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*B')]['perc'].sum(),
                'FB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('FB\*B')]['perc'].sum(),
            }
            Ratio_dic[_id] = _ratio

        _b = 'B3'
        _f = 'F4'
        Area_dic[_b] = sf[sf.MainGroup == _b].copy()
        Area_dic[_f] = sf[sf.MainGroup == _f].copy()
        _total = Area_dic[_b][area_used].sum() + Area_dic[_f][area_used].sum()
        Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
        Area_dic[_f]['perc'] = Area_dic[_f][area_used] / _total if _total != 0 else 0.0

        _ratio = {
            'F*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('F\*')]['perc'].sum(),
            'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
            'FF*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('FF\*')]['perc'].sum(),
            'FB*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('FB\*')]['perc'].sum(),
            'BF*': Area_dic[_f][Area_dic[_f]['Name'].str.contains('BF\*')]['perc'].sum(),
            'BB*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*')]['perc'].sum(),
            'FF*F': Area_dic[_f][Area_dic[_f]['Name'].str.contains('FF\*F')]['perc'].sum(),
            'BF*F': Area_dic[_f][Area_dic[_f]['Name'].str.contains('BF\*F')]['perc'].sum(),
            'FF*B': Area_dic[_f][Area_dic[_f]['Name'].str.contains('FF\*B')]['perc'].sum(),
            'BF*B': Area_dic[_f][Area_dic[_f]['Name'].str.contains('BF\*B')]['perc'].sum(),
            'BB*F': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*F')]['perc'].sum(),
            'FB*F': Area_dic[_b][Area_dic[_b]['Name'].str.contains('FB\*F')]['perc'].sum(),
            'BB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*B')]['perc'].sum(),
            'FB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('FB\*B')]['perc'].sum(),
        }
        Ratio_dic[3] = _ratio

        Ratio_dic['avg'] = {}
        for _key in keys:
            Ratio_dic['avg'][_key] = float(np.mean([Ratio_dic[1][_key], Ratio_dic[2][_key], Ratio_dic[3][_key]]))
        for _key in ['BB*F', 'FB*F', 'BB*B', 'FB*B']:
            Ratio_dic['avg'][_key] = float(np.mean([Ratio_dic[1][_key], Ratio_dic[2][_key]]))
        for _key in ['FF*F', 'BF*F', 'FF*B', 'BF*B']:
            Ratio_dic['avg'][_key] = float(np.mean([Ratio_dic[2][_key]]))

    elif profile.name == "PHBHHx":
        df['MainGroup'] = pd.Categorical(df['MainGroup'], ['B1X1', 'X3', 'B3', 'B2', 'X2', 'X4', 'B4', 'X5', 'X6'])
        sf = df[['MainGroup', 'Name', area_used]]
        Area_dic: dict[str, pd.DataFrame] = {}
        Area_dic['b1x1'] = sf[sf.MainGroup == 'B1X1'].copy()
        sum_b1x1 = Area_dic['b1x1'][area_used].sum()
        Area_dic['b1x1']['perc'] = Area_dic['b1x1'][area_used] / sum_b1x1 if sum_b1x1 != 0 else 0.0

        keys = ['X*', 'B*', 'X*X', 'X*B', 'B*X', 'B*B', 'XX*X', 'BX*X', 'XX*B', 'BX*B', 'BB*B', 'XB*B', 'BB*X', 'XB*X']
        str_keys = ['X\*', 'B\*', 'X\*X', 'X\*B', 'B\*X', 'B\*B', 'XX\*X', 'BX\*X', 'XX\*B', 'BX\*B', 'BB\*B', 'XB\*B', 'BB\*X', 'XB\*X']

        ratio = {}
        for _idx, _key in enumerate(keys):
            _value = np.nan
            if _key in ('X*X', 'B*X', 'X*B', 'B*B', 'B*', 'X*'):
                _value = Area_dic['b1x1'][Area_dic['b1x1']['Name'].str.contains(str_keys[_idx])]['perc'].sum()
            ratio[_key] = _value
        Ratio_dic[1] = ratio

        for _id in range(2, 3):
            _b = 'B' + str(_id)
            _x = 'X' + str(_id)
            Area_dic[_b] = sf[sf.MainGroup == _b].copy()
            Area_dic[_x] = sf[sf.MainGroup == _x].copy()
            _total = Area_dic[_b][area_used].sum() + Area_dic[_x][area_used].sum()
            Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
            Area_dic[_x]['perc'] = Area_dic[_x][area_used] / _total if _total != 0 else 0.0
            _ratio = {
                'X*': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*')]['perc'].sum(),
                'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
                'X*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*X')]['perc'].sum(),
                'X*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*B')]['perc'].sum(),
                'B*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*X')]['perc'].sum(),
                'B*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*B')]['perc'].sum(),
                'XX*X': np.nan, 'BX*X': np.nan, 'XX*B': np.nan, 'BX*B': np.nan,
                'BB*B': np.nan, 'XB*B': np.nan, 'BB*X': np.nan, 'XB*X': np.nan,
            }
            Ratio_dic[_id] = _ratio

        for _id in range(3, 4):
            _b = 'B' + str(_id)
            _x = 'X' + str(_id)
            Area_dic[_b] = sf[sf.MainGroup == _b].copy()
            Area_dic[_x] = sf[sf.MainGroup == _x].copy()
            _total = Area_dic[_b][area_used].sum() + Area_dic[_x][area_used].sum()
            Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
            Area_dic[_x]['perc'] = Area_dic[_x][area_used] / _total if _total != 0 else 0.0
            _ratio = {
                'X*': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*')]['perc'].sum(),
                'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
                'X*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*X')]['perc'].sum(),
                'X*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*B')]['perc'].sum(),
                'B*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*X')]['perc'].sum(),
                'B*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*B')]['perc'].sum(),
                'XX*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('XX\*X')]['perc'].sum(),
                'BX*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('BX\*X')]['perc'].sum(),
                'XX*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('XX\*B')]['perc'].sum(),
                'BX*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('BX\*B')]['perc'].sum(),
                'BB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*B')]['perc'].sum(),
                'XB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('XB\*B')]['perc'].sum(),
                'BB*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*X')]['perc'].sum(),
                'XB*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('XB\*X')]['perc'].sum(),
            }
            Ratio_dic[_id] = _ratio

        for _id in range(4, 5):
            _b = 'B' + str(_id)
            _x = 'X4'
            Area_dic[_b] = sf[sf.MainGroup == _b].copy()
            Area_dic[_x] = sf[sf.MainGroup == _x].copy()
            _total = Area_dic[_b][area_used].sum() + Area_dic[_x][area_used].sum()
            Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
            Area_dic[_x]['perc'] = Area_dic[_x][area_used] / _total if _total != 0 else 0.0
            _ratio = {
                'X*': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*')]['perc'].sum(),
                'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
                'X*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*X')]['perc'].sum(),
                'X*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*B')]['perc'].sum(),
                'B*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*X')]['perc'].sum(),
                'B*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*B')]['perc'].sum(),
                'XX*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('XX\*X')]['perc'].sum(),
                'BX*X': Area_dic[_x][Area_dic[_x]['Name'].str.contains('BX\*X')]['perc'].sum(),
                'XX*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('XX\*B')]['perc'].sum(),
                'BX*B': Area_dic[_x][Area_dic[_x]['Name'].str.contains('BX\*B')]['perc'].sum(),
                'BB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*B')]['perc'].sum(),
                'XB*B': Area_dic[_b][Area_dic[_b]['Name'].str.contains('XB\*B')]['perc'].sum(),
                'BB*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('BB\*X')]['perc'].sum(),
                'XB*X': Area_dic[_b][Area_dic[_b]['Name'].str.contains('XB\*X')]['perc'].sum(),
            }
            Ratio_dic[_id] = _ratio

        Ratio_dic['avg'] = {}
        for _key in keys:
            Ratio_dic['avg'][_key] = float(np.mean([
                Ratio_dic[1][_key], Ratio_dic[2][_key], Ratio_dic[3][_key], Ratio_dic[4][_key]
            ]))
        for _key in ['XX*X', 'BX*X', 'XX*B', 'BX*B']:
            Ratio_dic['avg'][_key] = float(np.mean(Ratio_dic[4][_key]))
        for _key in ['BB*B', 'XB*B', 'BB*X', 'XB*X']:
            Ratio_dic['avg'][_key] = float(np.mean([Ratio_dic[3][_key], Ratio_dic[4][_key]]))

    # Fit Kamiya statistical models (Models 1, 2, 3)
    avg_ratio = Ratio_dic['avg']
    Ratio_dic['model_1'] = compute_kamiya_model_1(avg_ratio, profile.keys_map)
    Ratio_dic['model_2'] = compute_kamiya_model_2(avg_ratio, profile.keys_map)
    Ratio_dic['model_3'] = compute_kamiya_model_3(avg_ratio, profile.keys_map)
    apply_copolymer_statistics(Ratio_dic, profile.keys_map)

    pha_df['Ratio_dic'] = Ratio_dic
    pha_df['peaks_data_df'] = df
    return pha_df


# =============================================================================
# MCL Terpolymer Ratios
# =============================================================================

def compute_area_and_ratio_monomers(pha_df: dict[str, Any], use_line_fitting: bool = False) -> dict[str, Any]:
    """Compute monomer proportions for PHBHHxHO_MCL."""
    peaks_data = pha_df.get('peaks_data', [])
    area_used = 'Area_curve' if use_line_fitting else 'Area_raw'

    if not peaks_data:
        pha_df['peaks_data_df'] = pd.DataFrame(
            columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve']
        )
        pha_df['Ratio_dic'] = {'avg': {'B*': 0, 'X*': 0, 'O': 0, 'X*X': 0, 'X*B': 0, 'X*O': 0, 'B*B': 0, 'B*X': 0, 'B*O': 0, 'O*O': 0, 'O*X': 0, 'O*B': 0}}
        return pha_df

    df = pd.DataFrame(peaks_data, columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve'])
    df['MainGroup'] = pd.Categorical(df['MainGroup'], [
        'B1X1O1', 'O3', 'X3', 'B3', 'B2', 'O2X2', 'X4', 'O4', 'O5', 'O6', 'O7', 'B4', 'X5', 'O8', 'X6'
    ])
    sf = df[['MainGroup', 'Name', area_used]]

    Area_dic: dict[str, pd.DataFrame] = {}
    Ratio_dic: dict[Any, dict[str, float]] = {}

    Area_dic['b1x1o1'] = sf[sf.MainGroup == 'B1X1O1'].copy()
    sum_b1 = Area_dic['b1x1o1'][area_used].sum()
    Area_dic['b1x1o1']['perc'] = Area_dic['b1x1o1'][area_used] / sum_b1 if sum_b1 != 0 else 0.0

    keys = ['X*', 'B*', 'O*']
    str_keys = ['X\*', 'B\*', 'O\*']
    ratio: dict[str, float] = {}
    for _idx, _key in enumerate(keys):
        ratio[_key] = Area_dic['b1x1o1'][Area_dic['b1x1o1']['Name'].str.contains(str_keys[_idx])]['perc'].sum()
    Ratio_dic[1] = ratio

    for _id in range(2, 3):
        _b = 'B' + str(_id)
        _xo = 'O2X2'
        Area_dic[_b] = sf[sf.MainGroup == _b].copy()
        Area_dic[_xo] = sf[sf.MainGroup == _xo].copy()
        _total = Area_dic[_b][area_used].sum() + Area_dic[_xo][area_used].sum()
        Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
        Area_dic[_xo]['perc'] = Area_dic[_xo][area_used] / _total if _total != 0 else 0.0
        Ratio_dic[_id] = {
            'X*': Area_dic[_xo][Area_dic[_xo]['Name'].str.contains('X\*')]['perc'].sum(),
            'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
            'O*': Area_dic[_xo][Area_dic[_xo]['Name'].str.contains('O\*')]['perc'].sum(),
        }

    for _id in range(3, 4):
        _b = 'B' + str(_id)
        _x = 'X' + str(_id)
        _o = 'O' + str(_id)
        Area_dic[_b] = sf[sf.MainGroup == _b].copy()
        Area_dic[_x] = sf[sf.MainGroup == _x].copy()
        Area_dic[_o] = sf[sf.MainGroup == _o].copy()
        _total = Area_dic[_b][area_used].sum() + Area_dic[_x][area_used].sum() + Area_dic[_o][area_used].sum()
        Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
        Area_dic[_x]['perc'] = Area_dic[_x][area_used] / _total if _total != 0 else 0.0
        Area_dic[_o]['perc'] = Area_dic[_o][area_used] / _total if _total != 0 else 0.0
        Ratio_dic[_id] = {
            'X*': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*')]['perc'].sum(),
            'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
            'O*': Area_dic[_o][Area_dic[_o]['Name'].str.contains('O\*')]['perc'].sum(),
        }

    for _id in range(4, 5):
        _b = 'B' + str(_id)
        _x = 'X' + str(_id)
        _o = 'O' + str(_id)
        Area_dic[_b] = sf[sf.MainGroup == _b].copy()
        Area_dic[_x] = sf[sf.MainGroup == _x].copy()
        Area_dic[_o] = sf[sf.MainGroup == _o].copy()
        _total = Area_dic[_b][area_used].sum() + Area_dic[_x][area_used].sum() + Area_dic[_o][area_used].sum()
        Area_dic[_b]['perc'] = Area_dic[_b][area_used] / _total if _total != 0 else 0.0
        Area_dic[_x]['perc'] = Area_dic[_x][area_used] / _total if _total != 0 else 0.0
        Area_dic[_o]['perc'] = Area_dic[_o][area_used] / _total if _total != 0 else 0.0
        Ratio_dic[_id] = {
            'X*': Area_dic[_x][Area_dic[_x]['Name'].str.contains('X\*')]['perc'].sum(),
            'B*': Area_dic[_b][Area_dic[_b]['Name'].str.contains('B\*')]['perc'].sum(),
            'O*': Area_dic[_o][Area_dic[_o]['Name'].str.contains('O\*')]['perc'].sum(),
        }

    Ratio_dic['avg'] = {}
    for _key in keys:
        Ratio_dic['avg'][_key] = float(np.mean([
            Ratio_dic[1][_key], Ratio_dic[2][_key], Ratio_dic[3][_key], Ratio_dic[4][_key]
        ]))

    pha_df['Ratio_dic'] = Ratio_dic
    pha_df['peaks_data_df'] = df
    return pha_df


def compute_area_and_ratio_monomers_mcl(pha_df: dict[str, Any], use_line_fitting: bool = False) -> dict[str, Any]:
    """Compute monomer proportions including mcl side chains for PHBHHxHO_MCL."""
    peaks_data = pha_df.get('peaks_data', [])
    area_used = 'Area_curve' if use_line_fitting else 'Area_raw'

    if not peaks_data:
        pha_df['peaks_data_df'] = pd.DataFrame(
            columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve']
        )
        pha_df['Ratio_dic_mcl'] = {'B*': 0, 'X*': 0, 'O': 0, 'X*X': 0, 'X*B': 0, 'X*O': 0, 'B*B': 0, 'B*X': 0, 'B*O': 0, 'O*O': 0, 'O*X': 0, 'O*B': 0, 'mcl*': 0}
        return pha_df

    df = pd.DataFrame(peaks_data, columns=['MainGroup', 'Name', 'Position', 'Height', 'Area_raw', 'Area_curve'])
    df['MainGroup'] = pd.Categorical(df['MainGroup'], [
        'B1X1O1', 'O3', 'X3', 'B3', 'B2', 'O2X2', 'X4', 'O4', 'O5', 'O6', 'O7', 'B4', 'X5', 'mcl', 'O8', 'X6'
    ])
    sf = df[['MainGroup', 'Name', area_used]]

    Area_dic: dict[str, pd.DataFrame] = {
        'B4': sf[sf.MainGroup == 'B4'].copy(),
        'mcl': sf[sf.MainGroup == 'mcl'].copy(),
        'O8': sf[sf.MainGroup == 'O8'].copy(),
        'X6': sf[sf.MainGroup == 'X6'].copy(),
    }
    _total = (
        Area_dic['B4'][area_used].sum()
        + Area_dic['mcl'][area_used].sum()
        + Area_dic['O8'][area_used].sum()
        + Area_dic['X6'][area_used].sum()
    )

    for k in ('B4', 'mcl', 'O8', 'X6'):
        Area_dic[k]['perc'] = Area_dic[k][area_used] / _total if _total != 0 else 0.0

    ratio = {
        'B*': Area_dic['B4'][Area_dic['B4']['Name'].str.contains('B\*')]['perc'].sum(),
        'mcl*': Area_dic['mcl'][Area_dic['mcl']['Name'].str.contains('mcl\*')]['perc'].sum(),
        'X*': Area_dic['X6'][Area_dic['X6']['Name'].str.contains('X\*')]['perc'].sum(),
        'O*': Area_dic['O8'][Area_dic['O8']['Name'].str.contains('O\*')]['perc'].sum(),
    }

    pha_df['Ratio_dic_mcl'] = ratio
    pha_df['peaks_data_df'] = df
    return pha_df


# =============================================================================
# Graph Generation
# =============================================================================

def proportion_graph_copolymer(
    x: list[float],
    y: list[float],
    pha_df: dict[str, Any],
    profile: PolymerProfile,
) -> tuple[Figure, Figure, Figure]:
    """Create proportion summary graphs: full spectrum, stacked bars, and carbon contributions."""
    Ratio_dic = pha_df['Ratio_dic']
    m_a, m_b = profile.monomers[1], profile.monomers[0]  # e.g. V, B or F, B or X, B
    km = profile.keys_map

    # Figure 1: Full spectrum
    fig = plt.figure(figsize=(8.27, 2))
    ax2 = fig.add_subplot()
    ax2.plot(x, y, linewidth=1)
    ax2.set_xlim((-5, 200))
    ax2.set_ylim((-0.1e12, 0.5e12))
    ax2.invert_xaxis()
    ax2.set_xlabel("ppm", fontsize=7)
    ax2.set_ylabel("Intensity", fontsize=7)

    # Figure 2: Stacked bar plot of monomers and dyads
    fig1 = plt.figure(figsize=(8.27, 3))
    ax = fig1.add_subplot()

    def plot_txt(rect, value, txt, xoffset, yoffset):
        if value > 0:
            rect = rect[0]
            yloc = rect.get_y() + rect.get_height() / 2.0 + 0.11 + yoffset * 0.18
            xloc = int(rect.get_width()) - 8 + xoffset * 3
            ax.text(xloc, yloc, '{0:0.1f}%'.format(value))
            yloc = rect.get_y() + rect.get_height() / 2.0 - 0.05 + yoffset * 0.18
            ax.text(xloc, yloc, txt)

    colors_g = sns.color_palette("Greens_r", 8)
    colors_b = sns.color_palette("Blues_r", 8)
    colors_r = sns.color_palette("YlOrRd", 5)

    _b = Ratio_dic['avg'][km['B*']] * 100
    _a = Ratio_dic['avg'][km['A*']] * 100
    _bb = Ratio_dic['avg'][km['BB']] * 100
    _ba = Ratio_dic['avg'][km['BA']] * 100
    _ab = Ratio_dic['avg'][km['AB']] * 100
    _aa = Ratio_dic['avg'][km['AA']] * 100

    offset = 1 if (_a < 5.5 or _b < 5.5) else 0
    rec = ax.barh(0, _a + _b, color=colors_g[2], align='center')
    plot_txt(rec, _a, km['A*'], 0, offset)
    rec = ax.barh(0, _b, color=colors_b[2], align='center')
    plot_txt(rec, _b, km['B*'], 0, -offset)

    offset = 1 if (_bb < 5.5 or _ba < 5.5 or _aa < 5.5 or _ab < 5.5) else 0
    rec = ax.barh(1, _bb + _ba + _ab + _aa, color=colors_g[4], align='center')
    plot_txt(rec, _aa, km['AA'], offset, -offset)
    rec = ax.barh(1, _bb + _ba + _ab, color=colors_g[3], align='center')
    plot_txt(rec, _ab, km['AB'], offset, offset)
    rec = ax.barh(1, _bb + _ba, color=colors_b[4], align='center')
    plot_txt(rec, _ba, km['BA'], -offset, offset)
    rec = ax.barh(1, _bb, color=colors_b[3], align='center')
    plot_txt(rec, _bb, km['BB'], -offset, -offset)

    # Figure 3: Carbon contribution bar chart
    fig2 = plt.figure(figsize=[8.27, 3])
    ax4 = fig2.add_subplot()

    df_dict = pd.DataFrame.from_dict(Ratio_dic).transpose().unstack()
    d = pd.DataFrame(df_dict).reset_index()
    d.columns = ['Resonances', 'Number', 'Pecentages']

    cat_order = [
        km['A*'], km['B*'], km['AA'], km['AB'], km['BA'], km['BB'],
        km['AAA'], km['BAA'], km['AAB'], km['BAB']
    ]
    d['Resonances'] = pd.Categorical(d['Resonances'], cat_order)
    d = d[d.Resonances != 'D']
    d = d[d.Number != 'model_1']
    d = d[d.Number != 'model_2']
    d = d[d.Number != 'model_3']

    ax4 = sns.barplot(x='Resonances', y='Pecentages', hue='Number', data=d, palette=colors_r)
    ax4.set_ylabel('')
    ax4.set_xlabel('')
    plt.legend(bbox_to_anchor=(0.0, 1.0, 1.0, 0.10), loc=3, ncol=8, mode="expand")
    plt.tight_layout()

    return fig, fig1, fig2


def proportion_graphs_mcl(
    x: list[float],
    y: list[float],
    pha_df_ratio: dict[str, Any],
    pha_df_ratio_mcl: dict[str, Any],
) -> tuple[Figure, Figure, Figure, Figure]:
    """Create summary graphs for PHBHHxHO_MCL."""
    Ratio_dic = pha_df_ratio['Ratio_dic']
    Ratio_dic_mcl = pha_df_ratio_mcl['Ratio_dic_mcl']

    fig = plt.figure(figsize=(8.27, 2))
    ax2 = fig.add_subplot()
    ax2.plot(x, y, linewidth=1)
    ax2.set_xlim((-5, 200))
    ax2.set_ylim((-0.1e12, 0.5e12))
    ax2.invert_xaxis()
    ax2.set_xlabel("ppm", fontsize=7)
    ax2.set_ylabel("Intensity", fontsize=7)

    fig1 = plt.figure(figsize=(8.27, 3))
    ax = fig1.add_subplot()

    colors_g = sns.color_palette("Greens_r", 8)
    colors_b = sns.color_palette("Blues_r", 8)
    colors_o = sns.color_palette("Oranges_r", 8)
    colors_r = sns.color_palette("YlOrRd", 5)

    _b = Ratio_dic['avg']['B*'] * 100
    _x = Ratio_dic['avg']['X*'] * 100
    _o = Ratio_dic['avg']['O*'] * 100

    ax.barh(0, _x + _b + _o, color=colors_o[2], align='center')
    ax.barh(0, _b + _x, color=colors_g[2], align='center')
    ax.barh(0, _b, color=colors_b[2], align='center')

    fig2 = plt.figure(figsize=[8.27, 3])
    ax4 = fig2.add_subplot()
    df_dict = pd.DataFrame.from_dict(Ratio_dic).transpose().unstack()
    d = pd.DataFrame(df_dict).reset_index()
    d.columns = ['Resonances', 'Number', 'Pecentages']
    d['Resonances'] = pd.Categorical(d['Resonances'], ['X*', 'B*', 'O*'])
    ax4 = sns.barplot(x='Resonances', y='Pecentages', hue='Number', data=d, palette=colors_r)
    ax4.set_ylabel('')
    ax4.set_xlabel('')
    plt.legend(bbox_to_anchor=(0.0, 1.0, 1.0, 0.10), loc=3, ncol=8, mode="expand")
    plt.tight_layout()

    fig3 = plt.figure(figsize=(8.27, 3))
    ax5 = fig3.add_subplot()
    df_dict_mcl = pd.DataFrame.from_dict(Ratio_dic_mcl, orient='index', columns=['Pecentages']).reset_index()
    df_dict_mcl.columns = ['Resonances', 'Pecentages']
    df_dict_mcl['Resonances'] = pd.Categorical(df_dict_mcl['Resonances'], ['B*', 'X*', 'O*', 'mcl*'])
    sns.barplot(x='Resonances', y='Pecentages', data=df_dict_mcl, palette=colors_r, ax=ax5)
    ax5.set_ylabel('')
    ax5.set_xlabel('')
    plt.tight_layout()

    return fig, fig1, fig2, fig3


# =============================================================================
# CSV Export
# =============================================================================

def export_csv_results(
    pha_df: dict[str, Any],
    file_path: str,
    profile: PolymerProfile,
    multiple_noise: Union[int, float] = 10,
) -> str:
    """Export calculated copolymer statistics and models to CSV."""
    Ratio_dic = pha_df['Ratio_dic']
    date = datetime.datetime.now()
    output_prefix = resolve_output_prefix(file_path)
    output_name = f"{output_prefix}_13C_n{multiple_noise}_Software_results_{date.strftime('%d-%m-%Y %H-%M-%S')}.csv"
    output_path = os.path.join(os.path.dirname(file_path), output_name)

    km = profile.keys_map

    with open(output_path, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow([output_name])
        spamwriter.writerow('')
        spamwriter.writerow(profile.csv_header)

        # Experimental rows (1, 2, 3, 4, avg)
        for _id in [1, 2, 3, 4, 'avg']:
            if _id not in Ratio_dic:
                continue
            r = Ratio_dic[_id]
            if profile.name == "PHBV":
                row = [
                    _id, r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_V'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*V'), r.get('V*'), r.get('V*B'), r.get('V*V'), '',
                    r.get('VV*V'), r.get('BV*V'), r.get('VV*B'), r.get('BV*B')
                ]
            elif profile.name == "P34HB":
                row = [
                    _id, r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_F'), '',
                    r.get('B*'), r.get('BB*'), r.get('FB*'), r.get('F*'), r.get('BF*'), r.get('FF*'), '',
                    r.get('FF*F'), r.get('BF*F'), r.get('FF*B'), r.get('BF*B'),
                    r.get('FB*F'), r.get('BB*F'), r.get('FB*B'), r.get('BB*B')
                ]
            elif profile.name == "PHBHHx":
                row = [
                    _id, r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_X'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*X'), r.get('X*'), r.get('X*B'), r.get('X*X'), '',
                    r.get('XX*X'), r.get('BX*X'), r.get('XX*B'), r.get('BX*B'),
                    r.get('XB*X'), r.get('BB*X'), r.get('XB*B'), r.get('BB*B')
                ]
            spamwriter.writerow(row)

        spamwriter.writerow('')

        # Model 1
        if 'model_1' in Ratio_dic:
            r = Ratio_dic['model_1']
            if profile.name == "PHBV":
                spamwriter.writerow([
                    'model_1', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_V'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*V'), r.get('V*'), r.get('V*B'), r.get('V*V'), '',
                    r.get('VV*V'), r.get('BV*V'), r.get('VV*B'), r.get('BV*B'),
                    r.get('VB*V'), r.get('BB*V'), r.get('VB*B'), r.get('BB*B')
                ])
            elif profile.name == "P34HB":
                spamwriter.writerow([
                    'model_1', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_F'), '',
                    r.get('B*'), r.get('BB*'), r.get('FB*'), r.get('F*'), r.get('BF*'), r.get('FF*'), '',
                    r.get('FF*F'), r.get('BF*F'), r.get('FF*B'), r.get('BF*B'),
                    r.get('FB*F'), r.get('BB*F'), r.get('FB*B'), r.get('BB*B')
                ])
            elif profile.name == "PHBHHx":
                spamwriter.writerow([
                    'model_1', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_X'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*X'), r.get('X*'), r.get('X*B'), r.get('X*X'), '',
                    r.get('XX*X'), r.get('BX*X'), r.get('XX*B'), r.get('BX*B'),
                    r.get('XB*X'), r.get('BB*X'), r.get('XB*B'), r.get('BB*B')
                ])

        # Model 2
        if 'model_2' in Ratio_dic:
            r = Ratio_dic['model_2']
            if profile.name == "PHBV":
                spamwriter.writerow([
                    'model_2', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_V'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*V'), r.get('V*'), r.get('V*B'), r.get('V*V'), '',
                    r.get('VV*V'), r.get('BV*V'), r.get('VV*B'), r.get('BV*B'),
                    r.get('VB*V'), r.get('BB*V'), r.get('VB*B'), r.get('BB*B'), '',
                    r.get('Pvv'), r.get('Pbb'), r.get('Pvb'), r.get('Pbv')
                ])
            elif profile.name == "P34HB":
                spamwriter.writerow([
                    'model_2', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_F'), '',
                    r.get('B*'), r.get('BB*'), r.get('FB*'), r.get('F*'), r.get('BF*'), r.get('FF*'), '',
                    r.get('FF*F'), r.get('BF*F'), r.get('FF*B'), r.get('BF*B'),
                    r.get('FB*F'), r.get('BB*F'), r.get('FB*B'), r.get('BB*B'), '',
                    r.get('Pff'), r.get('Pbb'), r.get('Pfb'), r.get('Pbf')
                ])
            elif profile.name == "PHBHHx":
                spamwriter.writerow([
                    'model_2', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_X'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*X'), r.get('X*'), r.get('X*B'), r.get('X*X'), '',
                    r.get('XX*X'), r.get('BX*X'), r.get('XX*B'), r.get('BX*B'),
                    r.get('XB*X'), r.get('BB*X'), r.get('XB*B'), r.get('BB*B'), '',
                    r.get('Pxx'), r.get('Pbb'), r.get('Pxb'), r.get('Pbx')
                ])

        # Model 3
        if 'model_3' in Ratio_dic:
            r = Ratio_dic['model_3']
            if profile.name == "PHBV":
                spamwriter.writerow([
                    'model_3', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_V'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*V'), r.get('V*'), r.get('V*B'), r.get('V*V'), '',
                    r.get('VV*V'), r.get('BV*V'), r.get('VV*B'), r.get('BV*B'),
                    r.get('VB*V'), r.get('BB*V'), r.get('VB*B'), r.get('BB*B'),
                    '', '', '', '', '', '', r.get('A'), r.get('B'), r.get('X')
                ])
            elif profile.name == "P34HB":
                spamwriter.writerow([
                    'model_3', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_F'), '',
                    r.get('B*'), r.get('BB*'), r.get('FB*'), r.get('F*'), r.get('BF*'), r.get('FF*'), '',
                    r.get('FF*F'), r.get('BF*F'), r.get('FF*B'), r.get('BF*B'),
                    r.get('FB*F'), r.get('BB*F'), r.get('FB*B'), r.get('BB*B'),
                    '', '', '', '', '', '', r.get('A'), r.get('B'), r.get('X')
                ])
            elif profile.name == "PHBHHx":
                spamwriter.writerow([
                    'model_3', r.get('D'), r.get('chi'), r.get('len_B'), r.get('len_X'), '',
                    r.get('B*'), r.get('B*B'), r.get('B*X'), r.get('X*'), r.get('X*B'), r.get('X*X'), '',
                    r.get('XX*X'), r.get('BX*X'), r.get('XX*B'), r.get('BX*B'),
                    r.get('XB*X'), r.get('BB*X'), r.get('XB*B'), r.get('BB*B'),
                    '', '', '', '', '', '', r.get('A'), r.get('B'), r.get('X')
                ])

        spamwriter.writerow('')
        for note in profile.csv_notes:
            spamwriter.writerow([note])

    return output_path


# =============================================================================
# PDF Report Generation
# =============================================================================

def create_pdf_report(
    file_path: str,
    figures: list[Figure],
    multiple_noise: Union[int, float],
    pha_df: Optional[dict[str, Any]] = None,
    profile: Optional[PolymerProfile] = None,
) -> str:
    """Generate high-resolution 2-page PDF report with embedded figures and metric tables."""
    date = datetime.datetime.now()
    sample_name = resolve_sample_name(file_path)
    output_prefix = resolve_output_prefix(file_path)
    output_pdf_name = f"{output_prefix}_13C_n{multiple_noise}_Software_results_{date.strftime('%d-%m-%Y %H-%M-%S')}.pdf"
    pdf_path = os.path.join(os.path.dirname(file_path), output_pdf_name)

    c = canvas.Canvas(pdf_path, pagesize=letter)
    c.drawCentredString(
        297.72, 750,
        f'NMR results for {sample_name} from PHACS-NMR software {date.strftime("%d-%m-%Y %H-%M-%S")} \n Noise factor: {multiple_noise}'
    )

    if profile and profile.has_mcl:
        # MCL 5-figure PDF layout
        for fig in figures:
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=300)
            buf.seek(0)
            if fig == figures[0]:
                c.drawString(20, 690, 'Figure 1: Full spectrum')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 520)
            elif fig == figures[1]:
                c.drawString(20, 470, 'Figure 2: Zoom on ranges of interest and peaks identification')
                img = Image(buf, width=595, height=400)
                img.drawOn(c, 10, 50)
            elif fig == figures[4]:
                c.drawString(20, 695, 'Figure 3: monomer class composition regarding side-chain length: look for monomer >C10')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 535)
            elif fig == figures[2]:
                c.drawString(20, 500, 'Figure 4: Calculation of monomers proportions and dyads proportions')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 340)
            elif fig == figures[3]:
                c.drawString(20, 310, 'Figure 5: Contribution of each couple of similar carbon on the average proportions')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 140)

            if fig == figures[1]:
                c.showPage()
                c.drawCentredString(
                    297.72, 750,
                    f'NMR results for {sample_name} from PHACS-NMR software {date.strftime("%d-%m-%Y %H-%M-%S")}'
                )
    else:
        # Standard copolymer 4-figure PDF layout
        Ratio_dic = pha_df['Ratio_dic'] if (pha_df and 'Ratio_dic' in pha_df) else {}
        avg = Ratio_dic.get('avg', {})

        m_a = profile.monomers[1] if (profile and len(profile.monomers) > 1) else 'A'
        len_a_key = profile.keys_map.get('len_A', f'len_{m_a}') if profile else f'len_{m_a}'

        for fig in figures:
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=300)
            buf.seek(0)

            if fig == figures[0]:
                c.drawString(20, 690, 'Figure 1: Full spectrum')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 520)
            elif fig == figures[1]:
                c.drawString(20, 470, 'Figure 2: Zoom on ranges of interest and peaks identification')
                img = Image(buf, width=595, height=400)
                img.drawOn(c, 10, 50)
            elif fig == figures[2]:
                c.drawString(20, 695, 'Figure 3: Calculation of monomers proportions and dyads proportions')
                img = Image(buf, width=595, height=175)
                img.drawOn(c, 0, 500)
            elif fig == figures[3]:
                c.drawString(20, 470, 'Figure 4: Contribution of each couple of similar carbon on the average proportions')
                img = Image(buf, width=595, height=150)
                img.drawOn(c, 0, 300)

                d_val = round(avg.get('D', 0.0), 2)
                chi_val = round(avg.get('chi', 0.0), 2)
                len_a_val = round(avg.get(len_a_key, 0.0), 1) if not np.isnan(avg.get(len_a_key, 0.0)) else 0.0
                len_b_val = round(avg.get('len_B', 0.0), 1) if not np.isnan(avg.get('len_B', 0.0)) else 0.0

                c.setFont("Helvetica-Bold", 14)
                c.drawString(50, 270, f"D = {d_val}")
                c.setFont("Helvetica", 12)
                c.drawString(50, 250, 'D = 1: Completely random copolymer,')
                c.drawString(50, 230, 'D > 1: Block copolymer or mix of different random copolymer,')
                c.drawString(50, 210, 'D < 1: Copolymer with an alternative tendency,')
                c.setFont("Helvetica-Bold", 14)
                c.drawString(50, 190, f"Chi = {chi_val}")
                c.setFont("Helvetica", 12)
                c.drawString(50, 170, 'chi = 1: Completely random copolymer,')
                c.drawString(50, 150, 'chi > 1: Copolymer with an alternating tendency,')
                c.drawString(50, 130, 'chi = 2: Completely alternating copolymer,')
                c.drawString(50, 110, 'chi < 1: Copolymer with a blocking tendency,')
                c.drawString(50, 90, 'chi = 0: Completely block copolymer')
                c.setFont("Helvetica-Bold", 14)
                c.drawString(50, 70, f"N({m_a}) = {len_a_val}; N(B) = {len_b_val}")
                c.setFont("Helvetica", 12)
                c.drawString(50, 50, "Number-average length for monomer blocks")

            if fig == figures[1]:
                c.showPage()
                c.drawCentredString(
                    297.72, 750,
                    f'NMR results for {sample_name} from PHACS-NMR software {date.strftime("%d-%m-%Y %H-%M-%S")}'
                )

    c.save()
    print(f"PDF report created: {output_pdf_name}")
    return pdf_path


# =============================================================================
# Master Analysis Pipeline
# =============================================================================

def run_analysis(
    file_path: str,
    profile: PolymerProfile,
    multiple_noise: Union[int, float] = 10,
    col_x_index: int = 3,
    col_y_index: int = 1,
    col_x_name: str = 'ppm',
    col_y_name: str = 'Intensity',
    header: Any = None,
    show_plot: bool = False,
) -> list[Figure]:
    """
    Run complete end-to-end analysis on NMR spectrum text file using specified profile.
    """
    x, y, pha_df, fig1 = split_graphic_regions(
        file_path,
        profile.spectral_regions,
        multiple_noise,
        profile.limits_builder,
        col_x_index=col_x_index,
        col_y_index=col_y_index,
        col_x_name=col_x_name,
        col_y_name=col_y_name,
        header=header,
        figsize=profile.figsize_regions,
        noise_as_int=profile.noise_as_int,
        monomer_colors=profile.monomer_colors,
    )

    if profile.has_mcl:
        pha_df_ratio = compute_area_and_ratio_monomers(pha_df)
        pha_df_ratio_mcl = compute_area_and_ratio_monomers_mcl(pha_df)
        fig, fig2, fig3, fig4 = proportion_graphs_mcl(x, y, pha_df_ratio, pha_df_ratio_mcl)
        figures = [fig, fig1, fig2, fig3, fig4]

        print("Creation of the pdf file")
        create_pdf_report(file_path, figures, multiple_noise, profile=profile)
        print("pdf file created")
    else:
        pha_df_s = compute_area_and_ratio(pha_df, profile)
        fig, fig2, fig3 = proportion_graph_copolymer(x, y, pha_df_s, profile)
        figures = [fig, fig1, fig2, fig3]

        print("Creation of the pdf file")
        create_pdf_report(file_path, figures, multiple_noise, pha_df=pha_df_s, profile=profile)

        if profile.has_csv:
            print("Creation of the excel file")
            export_csv_results(pha_df, file_path, profile, multiple_noise=multiple_noise)
            print("Excel file created")

    if show_plot:
        plt.show()

    return figures
