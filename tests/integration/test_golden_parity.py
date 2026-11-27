import io
import os
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import phacs_nmr
from tests.conftest import BENCHMARK_DATA_DIR as BENCHMARK_DIR


def parse_software_csv(csv_path: Path):
    """Parse software results CSV into a structured DataFrame of numeric rows."""
    # Rows 3-12 contain: id, D, chi, len_B, len_V, ..., model_1, model_2, model_3
    # Read ignoring initial title lines and footer documentation
    df = pd.read_csv(
        csv_path,
        skiprows=2,
        nrows=9,  # rows 1, 2, 3, 4, avg, empty, model_1, model_2, model_3
        index_col=0,
    )
    # Drop completely empty columns/rows if any
    df = df.dropna(how="all")
    return df


class TestGoldenParity:
    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_phbv_sample_hph18_parity(self, tmp_path, monkeypatch):
        """
        Verify that analyzing benchmark sample 'SI_13Cq_PHA HPH18 25 PHBV_11.txt'
        with noise=10 yields numerical parity with the historical baseline CSV.
        """
        input_txt = BENCHMARK_DIR / "SI_13Cq_PHA HPH18 25 PHBV_11.txt"
        historical_csv = BENCHMARK_DIR / "_SI_13Cq_PHA HPH18 25 PHBV_11_Software_results_26-06-2026 14-32-49.csv"

        assert input_txt.exists(), f"Input spectrum not found: {input_txt}"
        assert historical_csv.exists(), f"Historical CSV not found: {historical_csv}"

        # Copy input file to tmp_path to avoid creating test output artifacts in the benchmark folder
        # Preserve directory depth (>= 5 components) so legacy split('\\')[-5] does not fail
        nested_tmp = tmp_path / "data" / "samples" / "PHBV" / "run"
        nested_tmp.mkdir(parents=True)
        test_input = nested_tmp / input_txt.name
        test_input.write_bytes(input_txt.read_bytes())

        # Prevent matplotlib GUI popup during tests
        import matplotlib
        matplotlib.use("Agg")

        # Run analysis through phacs_nmr
        figures = phacs_nmr.start_analysis(str(test_input), multiple_noise=10)
        assert len(figures) > 0, "start_analysis returned no figures"

        # Find the newly generated CSV in nested_tmp
        generated_csvs = list(nested_tmp.glob("*_Software_results_*.csv"))
        assert len(generated_csvs) == 1, f"Expected 1 CSV, found {generated_csvs}"

        # Parse historical baseline and new generated CSV
        df_golden = parse_software_csv(historical_csv)
        df_new = parse_software_csv(generated_csvs[0])

        # Compare common numeric columns
        common_cols = [c for c in df_golden.columns if c in df_new.columns and not c.startswith("Unnamed")]
        assert len(common_cols) >= 5, f"Not enough common columns: {common_cols}"

        for row_idx in ["1", "2", "3", "4", "avg", "model_1", "model_2"]:
            if row_idx in df_golden.index and row_idx in df_new.index:
                for col in ["D", "chi", "B*", "V*", "B*B", "V*V"]:
                    if col in common_cols:
                        expected = float(df_golden.loc[row_idx, col])
                        actual = float(df_new.loc[row_idx, col])
                        if np.isnan(expected):
                            assert np.isnan(actual), f"Row {row_idx} col {col}: expected NaN, got {actual}"
                        else:
                            assert np.isclose(actual, expected, rtol=1e-4, atol=1e-4), (
                                f"Row {row_idx} col {col}: expected {expected}, got {actual}"
                            )

    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_p34hb_sample_parity(self, tmp_path):
        """
        Verify that analyzing benchmark sample 'SI_13Cq_PHA A1000P 50mg_mL P34HB_12.txt'
        with noise=10 yields numerical parity with historical baseline CSV.
        """
        input_txt = BENCHMARK_DIR / "SI_13Cq_PHA A1000P 50mg_mL P34HB_12.txt"
        historical_csv = BENCHMARK_DIR / "_SI_13Cq_PHA A1000P 50mg_mL P34HB_12_Software_results_26-06-2026 14-32-26.csv"

        assert input_txt.exists() and historical_csv.exists()

        nested_tmp = tmp_path / "data" / "samples" / "P34HB" / "run"
        nested_tmp.mkdir(parents=True)
        test_input = nested_tmp / input_txt.name
        test_input.write_bytes(input_txt.read_bytes())

        import matplotlib
        matplotlib.use("Agg")

        figures = phacs_nmr.start_analysis(str(test_input), multiple_noise=10)
        assert len(figures) > 0

        generated_csvs = list(nested_tmp.glob("*_Software_results_*.csv"))
        assert len(generated_csvs) == 1

        df_golden = parse_software_csv(historical_csv)
        df_new = parse_software_csv(generated_csvs[0])

        common_cols = [c for c in df_golden.columns if c in df_new.columns and not c.startswith("Unnamed")]
        for row_idx in ["1", "2", "3", "4", "avg"]:
            if row_idx in df_golden.index and row_idx in df_new.index:
                for col in ["B*", "F*", "B*B", "B*F", "F*B"]:
                    if col in common_cols:
                        expected = float(df_golden.loc[row_idx, col])
                        actual = float(df_new.loc[row_idx, col])
                        if np.isnan(expected):
                            assert np.isnan(actual)
                        else:
                            assert np.isclose(actual, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_phbhhx_sample_parity(self, tmp_path):
        """
        Verify that analyzing benchmark sample 'SI_13Cq_PHA BP 330-05 PHBHHx_11.txt'
        with noise=10 yields numerical parity with historical baseline CSV.
        """
        input_txt = BENCHMARK_DIR / "SI_13Cq_PHA BP 330-05 PHBHHx_11.txt"
        historical_csv = BENCHMARK_DIR / "_SI_13Cq_PHA BP 330-05 PHBHHx_11_Software_results_26-06-2026 14-32-35.csv"

        assert input_txt.exists() and historical_csv.exists()

        nested_tmp = tmp_path / "data" / "samples" / "PHBHHx" / "run"
        nested_tmp.mkdir(parents=True)
        test_input = nested_tmp / input_txt.name
        test_input.write_bytes(input_txt.read_bytes())

        import matplotlib
        matplotlib.use("Agg")

        figures = phacs_nmr.start_analysis(str(test_input), multiple_noise=10)
        assert len(figures) > 0

        generated_csvs = list(nested_tmp.glob("*_Software_results_*.csv"))
        assert len(generated_csvs) == 1

        df_golden = parse_software_csv(historical_csv)
        df_new = parse_software_csv(generated_csvs[0])

        common_cols = [c for c in df_golden.columns if c in df_new.columns and not c.startswith("Unnamed")]
        for row_idx in ["1", "2", "3", "4", "avg"]:
            if row_idx in df_golden.index and row_idx in df_new.index:
                for col in ["D", "chi", "B*", "X*", "B*B", "X*X"]:
                    if col in common_cols:
                        expected = float(df_golden.loc[row_idx, col])
                        actual = float(df_new.loc[row_idx, col])
                        if np.isnan(expected):
                            assert np.isnan(actual)
                        else:
                            assert np.isclose(actual, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_homopolymer_early_exit(self, tmp_path):
        """
        Verify that analyzing pure PHB homopolymer 'SI_13Cq_PHA Y1000P PHBV_11.txt'
        correctly classifies as homopolymer [0,1,0,0,0,0] and cleanly returns empty figures.
        """
        input_txt = BENCHMARK_DIR / "SI_13Cq_PHA Y1000P PHBV_11.txt"
        assert input_txt.exists()

        nested_tmp = tmp_path / "data" / "samples" / "Homopolymer" / "run"
        nested_tmp.mkdir(parents=True)
        test_input = nested_tmp / input_txt.name
        test_input.write_bytes(input_txt.read_bytes())

        figures = phacs_nmr.start_analysis(str(test_input), multiple_noise=10)
        assert figures == [], "Expected empty figures list for PHB homopolymer"

    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_a1000p_11_robust_fitting(self, tmp_path):
        """
        Verify robust lineshape fitting on previously failing 'SI_13Cq_PHA A1000P 50mg_mL P34HB_11.txt'.
        Previously crashed with 'cannot convert float NaN to integer' due to inner loop variable shadowing.
        Now completes successfully and produces 4 figures and results.
        """
        input_txt = BENCHMARK_DIR / "SI_13Cq_PHA A1000P 50mg_mL P34HB_11.txt"
        assert input_txt.exists()

        nested_tmp = tmp_path / "data" / "samples" / "P34HB_recovered" / "run"
        nested_tmp.mkdir(parents=True)
        test_input = nested_tmp / input_txt.name
        test_input.write_bytes(input_txt.read_bytes())

        figures = phacs_nmr.start_analysis(str(test_input), multiple_noise=10)
        assert len(figures) == 4, f"Expected 4 figures, got {len(figures)}"

