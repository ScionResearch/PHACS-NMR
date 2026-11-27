"""
Unit tests for core.io file parsing, error robustness, and monomer sensitivity filtering.
Verifies read_ascii_spectrum, read_dcon_file, find_monomers for all 9 PHA monomers,
and boundary threshold behavior in results_dconv without relying on historical outputs.
"""
from pathlib import Path
import pytest
import numpy as np

from phacs_nmr.core.io import (
    read_ascii_spectrum,
    read_dcon_file,
    find_monomers,
    results_dconv,
)


def create_mock_dconpeaks_file(filepath: Path, peak_records: list[tuple[float, float, float]]):
    """
    Helper to create synthetic TopSpin dconpeaks.txt format.
    peak_records: list of (ppm, intensity, area)
    """
    lines = [
        "TopSpin Deconv Output\n",
        "Date\n",
        "Sample\n",
        "Param 1\n",
        "Param 2\n",
        "Param 3\n",
        "Param 4\n",
        "Param 5\n",
    ]
    for ppm, intensity, area in peak_records:
        lines.append(f"{ppm:.4f} 0 0 0 {intensity:.1f} {area:.1f} 0.01\n")
    filepath.parent.mkdir(parents=True, exist_ok=True)
    filepath.write_text("".join(lines), encoding="cp1252")


class TestCoreIoRobustness:
    """Test suite for I/O parsing robustness and boundary filtering."""

    # -------------------------------------------------------------------------
    # 1. read_ascii_spectrum error handling & malformed input
    # -------------------------------------------------------------------------

    def test_read_ascii_spectrum_non_existent_file(self, tmp_path):
        """Non-existent file must raise FileNotFoundError."""
        missing = tmp_path / "non_existent.txt"
        with pytest.raises(FileNotFoundError):
            read_ascii_spectrum(missing)

    def test_read_ascii_spectrum_skips_malformed_and_blank_lines(self, tmp_path):
        """Malformed lines (< 4 columns), comments, and empty lines must be safely ignored."""
        content = (
            "# TopSpin convbin2asc export\n"
            "\n"
            "Corrupt line with three cols 1 2\n"
            "1, 100.5, 0.0, 169.5\n"
            "2  200.5  0.0  169.4\n"
            "Another bad line\n"
            "3, 300.5, 0.0, 169.3\n"
        )
        spec_file = tmp_path / "noisy_ascii.txt"
        spec_file.write_text(content, encoding="cp1252")

        x, y = read_ascii_spectrum(spec_file)
        assert len(x) == 3
        assert len(y) == 3
        assert x == [169.5, 169.4, 169.3]
        assert y == [100.5, 200.5, 300.5]

    # -------------------------------------------------------------------------
    # 2. read_dcon_file header skipping & standard line filtering
    # -------------------------------------------------------------------------

    def test_read_dcon_file_short_file(self, tmp_path):
        """File with fewer than 8 header lines must return empty lists."""
        short_file = tmp_path / "short_dcon.txt"
        short_file.write_text("Line1\nLine2\nLine3\n", encoding="cp1252")

        result = read_dcon_file(short_file)
        assert result["peaks_ppm"] == []
        assert result["peaks_intensity"] == []
        assert result["peaks_area"] == []

    def test_read_dcon_file_skips_std_and_incomplete_rows(self, tmp_path):
        """Rows starting with 'STD' or having < 7 items must be ignored."""
        lines = [
            f"Header {i}\n" for i in range(8)
        ] + [
            "STD internal reference peak 1 2 3 4 5 6 7\n",
            "5.2620 0 0 0 1000 5000 0.01\n",
            "INCOMPLETE 1 2\n",
            "1.2700 0 0 0 2000 10000 0.01\n",
        ]
        dcon_file = tmp_path / "dcon_with_std.txt"
        dcon_file.write_text("".join(lines), encoding="cp1252")

        result = read_dcon_file(dcon_file)
        assert len(result["peaks_ppm"]) == 2
        assert float(result["peaks_ppm"][0]) == pytest.approx(5.2620)
        assert float(result["peaks_ppm"][1]) == pytest.approx(1.2700)

    # -------------------------------------------------------------------------
    # 3. find_monomers recognition across all 9 PHA monomers
    # -------------------------------------------------------------------------

    def test_find_monomers_all_mcl_monomers(self):
        """
        Verify that find_monomers correctly identifies all MCL monomers when 4HB is absent:
        [B, V, X, Hp, O, N, >C9, >C9_1] (p_id[1:] are all 1, p_id[0] == 0)
        """
        # Diagnostic chemical shifts from limits dictionary:
        # B: 5.2620 (Est), 1.2700 (Met)
        # V: 0.8950 (Met), 5.1530 (Est)
        # X: 5.2050 (Est), 0.9110 (Met)
        # Hp: 0.9010 (Met)
        # O: 5.1820 (Est), 0.8790 (Met)
        # N: 0.8760 (Met)
        # >C9: 0.8769 (Met)
        # >C9_1: 0.8815 (Met)
        all_ppms = [
            5.2620, 1.2700, 0.8950, 5.1530,
            5.2050, 0.9110, 0.9010, 5.1820, 0.8790,
            0.8760, 0.8769, 0.8815
        ]
        dconv_data = {
            "peaks_ppm": [str(p) for p in all_ppms],
            "peaks_intensity": [1000.0] * len(all_ppms),
            "peaks_area": [5000.0] * len(all_ppms),
        }

        p_id, est_area, est_int = find_monomers(dconv_data)

        # Expected: F=0 (4HB absent), and all 8 other monomers flagged as 1
        assert p_id[0] == 0  # F
        assert p_id[1:] == [1, 1, 1, 1, 1, 1, 1, 1]

    def test_find_monomers_scl_with_4hb(self):
        """
        Verify that find_monomers routes through the 4HB/CH2 pathway when 4HB is present:
        [F, B, V] (p_id[:3] == [1, 1, 1])
        """
        peaks = [
            4.1085,  # F Est
            2.3520,  # F CH2
            5.2620,  # B Est
            2.5975,  # B CH2
            5.1530,  # V Est
            2.5580,  # V CH2
        ]
        dconv_data = {
            "peaks_ppm": [str(p) for p in peaks],
            "peaks_intensity": [20000.0] * len(peaks),
            "peaks_area": [50000.0] * len(peaks),
        }

        p_id, est_area, est_int = find_monomers(dconv_data)
        assert p_id[:3] == [1, 1, 1]

    # -------------------------------------------------------------------------
    # 4. results_dconv sensitivity boundary filtering
    # -------------------------------------------------------------------------

    def test_sensitivity_boundary_filtering(self, tmp_path):
        """
        Verify monomer sensitivity filtering right at the boundary:
        - Peak with relative intensity 0.008 (0.8% < 1.0% threshold) is excluded.
        - Peak with relative intensity 0.015 (1.5% > 1.0% threshold) is included.
        """
        dcon_file = tmp_path / "pdata" / "1" / "dconpeaks.txt"

        # Case 1: Monomer V has intensity 800 vs B intensity 100,000 (0.8% < 1%)
        peaks_below = [
            (5.2620, 100000.0, 500000.0),  # B Est
            (1.2700, 100000.0, 500000.0),  # B Met
            (0.8950, 800.0, 4000.0),       # V Met (0.8%)
        ]
        create_mock_dconpeaks_file(dcon_file, peaks_below)
        headers, p_id_below, p_type_below, _, _ = results_dconv(str(dcon_file), sensitivity=0.01)

        # V (index 2) must be filtered out
        assert p_id_below[1] == 1  # B retained
        assert p_id_below[2] == 0  # V excluded (< 1%)
        assert p_type_below == "PH-B"

        # Case 2: Monomer V has intensity 1500 vs B intensity 100,000 (1.5% > 1%)
        peaks_above = [
            (5.2620, 100000.0, 500000.0),  # B Est
            (1.2700, 100000.0, 500000.0),  # B Met
            (0.8950, 1500.0, 7500.0),      # V Met (1.5%)
        ]
        create_mock_dconpeaks_file(dcon_file, peaks_above)
        headers, p_id_above, p_type_above, _, _ = results_dconv(str(dcon_file), sensitivity=0.01)

        # V (index 2) must be retained
        assert p_id_above[1] == 1  # B retained
        assert p_id_above[2] == 1  # V retained (> 1%)
        assert p_type_above == "PH-B-V"
