"""
Unit tests for core.engine ratio and area computations.
Tests compute_area_and_ratio across PHBV, P34HB, and PHBHHx profiles,
as well as MCL terpolymer ratio functions, using in-memory mock peak data.
"""
from pathlib import Path
import pytest
import numpy as np
import pandas as pd

from phacs_nmr.core.engine import (
    compute_area_and_ratio,
    compute_area_and_ratio_monomers,
    compute_area_and_ratio_monomers_mcl,
    export_csv_results,
)
from phacs_nmr.core.config import PHBV_PROFILE, P34HB_PROFILE, PHBHHX_PROFILE


class TestCoreEngineRatios:
    """Test suite for compute_area_and_ratio logic and edge cases."""

    def test_empty_peaks_data_handling(self):
        """Empty peaks_data must return default structures with zeros and not raise exceptions."""
        for profile in [PHBV_PROFILE, P34HB_PROFILE, PHBHHX_PROFILE]:
            pha_df = {"peaks_data": []}
            result = compute_area_and_ratio(pha_df, profile)

            assert "peaks_data_df" in result
            assert isinstance(result["peaks_data_df"], pd.DataFrame)
            assert len(result["peaks_data_df"]) == 0
            assert "Ratio_dic" in result
            assert "avg" in result["Ratio_dic"]
            assert result["Ratio_dic"]["avg"]["D"] == 0
            assert result["Ratio_dic"]["avg"]["chi"] == 0

    def test_phbv_synthetic_peaks_ratio_calculation(self):
        """
        Verify PHBV ratio calculation from synthetic peak records.
        Tests area integration across B1V1, V2/B2, V3/B3, V4/B4 regions,
        verifying row generation (1, 2, 3, 4, avg) and Kamiya models (model_1, model_2, model_3).
        """
        # Build synthetic peaks_data: [MainGroup, Name, Position, Height, Area_raw, Area_curve]
        # True proportions: ~70% B and ~30% V
        peaks_data = [
            # Region 1: B1V1
            ["B1V1", "B*B", 169.15, 1000.0, 4900.0, 4900.0],
            ["B1V1", "B*V", 169.35, 500.0, 2100.0, 2100.0],
            ["B1V1", "V*B", 169.31, 500.0, 2100.0, 2100.0],
            ["B1V1", "V*V", 169.50, 200.0, 900.0, 900.0],
            ["B1V1", "B*", 169.20, 1500.0, 7000.0, 7000.0],
            ["B1V1", "V*", 169.40, 700.0, 3000.0, 3000.0],
            # Region 2: B2 & V2
            ["B2", "B*", 40.80, 1000.0, 7000.0, 7000.0],
            ["B2", "B*B", 40.82, 700.0, 4900.0, 4900.0],
            ["B2", "B*V", 40.78, 300.0, 2100.0, 2100.0],
            ["V2", "V*", 38.70, 500.0, 3000.0, 3000.0],
            ["V2", "V*B", 38.80, 350.0, 2100.0, 2100.0],
            ["V2", "V*V", 38.65, 150.0, 900.0, 900.0],
            ["V2", "VV*V", 38.65, 50.0, 270.0, 270.0],
            ["V2", "BV*V", 38.69, 100.0, 630.0, 630.0],
            ["V2", "VV*B", 38.78, 100.0, 630.0, 630.0],
            ["V2", "BV*B", 38.82, 250.0, 1470.0, 1470.0],
            # Region 3: B3 & V3
            ["B3", "B*", 67.60, 1000.0, 7000.0, 7000.0],
            ["B3", "B*B", 67.60, 700.0, 4900.0, 4900.0],
            ["B3", "B*V", 67.68, 300.0, 2100.0, 2100.0],
            ["V3", "V*", 71.90, 500.0, 3000.0, 3000.0],
            ["V3", "V*B", 71.88, 350.0, 2100.0, 2100.0],
            ["V3", "V*V", 71.93, 150.0, 900.0, 900.0],
            # Region 4: B4 & V4
            ["B4", "B*", 19.78, 1000.0, 7000.0, 7000.0],
            ["B4", "B*B", 19.77, 700.0, 4900.0, 4900.0],
            ["B4", "B*V", 19.73, 300.0, 2100.0, 2100.0],
            ["V4", "V*", 26.80, 500.0, 3000.0, 3000.0],
            ["V4", "V*B", 26.86, 350.0, 2100.0, 2100.0],
            ["V4", "V*V", 26.76, 150.0, 900.0, 900.0],
            ["V4", "VV*V", 26.76, 50.0, 270.0, 270.0],
            ["V4", "BV*V", 26.80, 100.0, 630.0, 630.0],
            ["V4", "VV*B", 26.83, 100.0, 630.0, 630.0],
            ["V4", "BV*B", 26.88, 250.0, 1470.0, 1470.0],
        ]

        pha_df = {"peaks_data": peaks_data}
        result = compute_area_and_ratio(pha_df, PHBV_PROFILE, use_line_fitting=False)

        ratio_dic = result["Ratio_dic"]
        # Required rows
        for row in [1, 2, 3, 4, "avg", "model_1", "model_2", "model_3"]:
            assert row in ratio_dic, f"Missing row {row} in Ratio_dic"

        # Check average proportions ~ 0.70 B* and ~ 0.30 V*
        avg = ratio_dic["avg"]
        assert np.isclose(avg["B*"] + avg["V*"], 1.0, atol=0.05)
        assert 0.65 <= avg["B*"] <= 0.75
        assert 0.25 <= avg["V*"] <= 0.35

        # Check Kamiya D and chi for model_1 are 1.0
        assert np.isclose(ratio_dic["model_1"]["D"], 1.0, atol=1e-4)
        assert np.isclose(ratio_dic["model_1"]["chi"], 1.0, atol=1e-4)

    def test_p34hb_synthetic_peaks_ratio_calculation(self):
        """Verify P34HB ratio calculation for B and F monomers."""
        peaks_data = [
            # B1 & F1
            ["B1", "B*", 169.2, 1000.0, 6000.0, 6000.0],
            ["B1", "BB*B", 169.15, 600.0, 3600.0, 3600.0],
            ["B1", "FB*B", 169.25, 400.0, 2400.0, 2400.0],
            ["F1", "F*", 172.0, 700.0, 4000.0, 4000.0],
            ["F1", "F*B", 171.9, 700.0, 4000.0, 4000.0],
            # B2 & F2
            ["B2", "B*", 40.8, 1000.0, 6000.0, 6000.0],
            ["F2", "F*", 30.7, 700.0, 4000.0, 4000.0],
            # B3 & F4
            ["B3", "B*", 67.6, 1000.0, 6000.0, 6000.0],
            ["F4", "F*", 63.5, 700.0, 4000.0, 4000.0],
            # B4 & F3
            ["B4", "B*", 19.8, 1000.0, 6000.0, 6000.0],
            ["F3", "F*", 24.0, 700.0, 4000.0, 4000.0],
        ]

        pha_df = {"peaks_data": peaks_data}
        result = compute_area_and_ratio(pha_df, P34HB_PROFILE, use_line_fitting=False)
        ratio_dic = result["Ratio_dic"]

        assert "avg" in ratio_dic
        avg = ratio_dic["avg"]
        assert np.isclose(avg["B*"] + avg["F*"], 1.0, atol=0.05)
        assert 0.55 <= avg["B*"] <= 0.65
        assert 0.35 <= avg["F*"] <= 0.45

    def test_phbhhx_synthetic_peaks_ratio_calculation(self):
        """Verify PHBHHx ratio calculation for B and X monomers."""
        peaks_data = [
            # B1X1
            ["B1X1", "B*", 169.2, 1000.0, 8000.0, 8000.0],
            ["B1X1", "X*", 169.4, 300.0, 2000.0, 2000.0],
            ["B1X1", "B*B", 169.15, 600.0, 6400.0, 6400.0],
            ["B1X1", "B*X", 169.3, 200.0, 1600.0, 1600.0],
            ["B1X1", "X*B", 169.35, 200.0, 1600.0, 1600.0],
            ["B1X1", "X*X", 169.5, 100.0, 400.0, 400.0],
            # B2 & X2
            ["B2", "B*", 40.8, 1000.0, 8000.0, 8000.0],
            ["X2", "X*", 39.0, 300.0, 2000.0, 2000.0],
            # B3 & X3
            ["B3", "B*", 67.6, 1000.0, 8000.0, 8000.0],
            ["X3", "X*", 70.6, 300.0, 2000.0, 2000.0],
            # B4 & X4
            ["B4", "B*", 19.8, 1000.0, 8000.0, 8000.0],
            ["X4", "X*", 36.0, 300.0, 2000.0, 2000.0],
        ]

        pha_df = {"peaks_data": peaks_data}
        result = compute_area_and_ratio(pha_df, PHBHHX_PROFILE, use_line_fitting=False)
        ratio_dic = result["Ratio_dic"]

        assert "avg" in ratio_dic
        avg = ratio_dic["avg"]
        assert np.isclose(avg["B*"] + avg["X*"], 1.0, atol=0.05)
        assert 0.75 <= avg["B*"] <= 0.85
        assert 0.15 <= avg["X*"] <= 0.25

    def test_mcl_monomer_ratios_empty_and_synthetic(self):
        """Verify MCL monomer ratio functions for empty input and valid peaks."""
        # 1. Empty input
        empty_res = compute_area_and_ratio_monomers({"peaks_data": []})
        assert "Ratio_dic" in empty_res
        assert empty_res["Ratio_dic"]["avg"]["B*"] == 0

        empty_mcl = compute_area_and_ratio_monomers_mcl({"peaks_data": []})
        assert "Ratio_dic_mcl" in empty_mcl
        assert empty_mcl["Ratio_dic_mcl"]["B*"] == 0

        # 2. Synthetic terpolymer data (B, X, O)
        mcl_peaks = [
            ["B1X1O1", "B*", 169.2, 1000.0, 5000.0, 5000.0],
            ["B1X1O1", "X*", 169.4, 600.0, 3000.0, 3000.0],
            ["B1X1O1", "O*", 169.5, 400.0, 2000.0, 2000.0],
            ["B2", "B*", 40.8, 1000.0, 5000.0, 5000.0],
            ["O2X2", "O*", 39.1, 400.0, 2000.0, 2000.0],
            ["O2X2", "X*", 39.0, 600.0, 3000.0, 3000.0],
        ]
        res = compute_area_and_ratio_monomers({"peaks_data": mcl_peaks})
        assert "Ratio_dic" in res
        ratio_1 = res["Ratio_dic"][1]
        assert np.isclose(ratio_1["B*"] + ratio_1["X*"] + ratio_1["O*"], 1.0, atol=1e-3)
        assert np.isclose(ratio_1["B*"], 0.5, atol=1e-3)
        assert np.isclose(ratio_1["X*"], 0.3, atol=1e-3)
        assert np.isclose(ratio_1["O*"], 0.2, atol=1e-3)

    @pytest.mark.parametrize("profile", [PHBV_PROFILE, P34HB_PROFILE, PHBHHX_PROFILE])
    def test_export_csv_results(self, tmp_path, profile):
        """Verify export_csv_results generates well-structured CSV file with header, rows, models, and notes."""
        sample_path = tmp_path / "data" / "run" / "sample1" / "10" / "pdata" / "1" / "ascii-spec.txt"
        sample_path.parent.mkdir(parents=True, exist_ok=True)
        sample_path.write_text("1, 100, 0, 169.5\n", encoding="cp1252")

        # Fake Ratio_dic
        ratio_dic = {
            1: {"D": 1.0, "chi": 1.0, "B*": 0.7, "V*": 0.3, "F*": 0.3, "X*": 0.3},
            "avg": {"D": 1.0, "chi": 1.0, "B*": 0.7, "V*": 0.3, "F*": 0.3, "X*": 0.3},
            "model_1": {"D": 1.0, "chi": 1.0, "B*": 0.7, "V*": 0.3, "F*": 0.3, "X*": 0.3},
            "model_2": {"D": 1.0, "chi": 1.0, "B*": 0.7, "V*": 0.3, "F*": 0.3, "X*": 0.3},
            "model_3": {"D": 1.0, "chi": 1.0, "B*": 0.7, "V*": 0.3, "F*": 0.3, "X*": 0.3, "A": 0.7, "B": 0.2, "X": 0.6},
        }
        pha_df = {"Ratio_dic": ratio_dic}

        csv_path = export_csv_results(pha_df, sample_path, profile, multiple_noise=10)
        assert Path(csv_path).exists()
        content = Path(csv_path).read_text(encoding="utf-8")
        assert "Software_results" in csv_path
        assert "model_1" in content
        assert "model_2" in content
        assert "model_3" in content
        assert "Sources:" in content
