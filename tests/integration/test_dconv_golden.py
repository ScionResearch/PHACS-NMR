from pathlib import Path
import pytest
import numpy as np
import pandas as pd
from phacs_nmr.core import io as Read_dcon_file
from tests.conftest import GLOBAL_SOFTWARE_DIR


def parse_dconv_results_csv(csv_path: Path):
    """Parse dconpeaks software results CSV to extract monomer intensities and areas."""
    lines = csv_path.read_text(encoding="cp1252", errors="ignore").splitlines()
    data = {}
    for i, line in enumerate(lines):
        if line.startswith("Method,") or (len(line.split(",")) > 1 and "4-HydroxyButyrate" in line):
            header = [c.strip() for c in line.split(",")]
            data["headers"] = header[1:]
        elif line.startswith("Intensities,"):
            vals = [float(v.strip()) for v in line.split(",")[1:]]
            data["intensities"] = vals
        elif line.startswith("Area,"):
            vals = [float(v.strip()) for v in line.split(",")[1:]]
            data["areas"] = vals
    return data


class TestDconvGoldenParity:
    @pytest.mark.skipif(GLOBAL_SOFTWARE_DIR is None or not GLOBAL_SOFTWARE_DIR.exists(), reason="Global software benchmark directory not configured or found")
    def test_phbv_dconpeaks_golden_match(self, tmp_path):
        """
        Verify that analyzing real PHBV dconpeaks.txt (250818_CF_HPH18 25 PHBV)
        produces exact parity with historical 18-11-2025 results CSV.
        """
        sample_dir = GLOBAL_SOFTWARE_DIR / "250818_CF_HPH18 25 PHBV" / "10" / "pdata" / "1"
        dcon_file = sample_dir / "dconpeaks.txt"
        historical_csv = sample_dir / "250818_CF_HPH18 25 PHBV_dconpeaks_Software_results_18-11-2025 14-18-37.csv"

        assert dcon_file.exists(), f"File not found: {dcon_file}"
        assert historical_csv.exists(), f"File not found: {historical_csv}"

        # Copy to isolated tmp_path preserving depth >= 5 components
        nested_tmp = tmp_path / "global" / "sample" / "250818_CF_HPH18 25 PHBV" / "10" / "pdata" / "1"
        nested_tmp.mkdir(parents=True)
        test_file = nested_tmp / dcon_file.name
        test_file.write_bytes(dcon_file.read_bytes())

        # Sensitivity = 0.01 matches original analysis
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(test_file), sensitivity=0.01)

        # In baseline CSV:
        # PH-B-V
        # Method,4-HydroxyButyrate,3-HydroxyButyrate,3-HydroxyValerate,3-HydroxyHexanoate,3-HydroxyOctanoate
        # Intensities,0.0,0.796,0.204,0.0,0.0
        assert polymer_type == "PH-B-V"
        assert p_id[1] == 1  # 3HB
        assert p_id[2] == 1  # 3HV
        assert p_id[0] == 0  # 4HB

        # In the 2026 model, the 3HV fraction is scaled by 0.8401:
        # Pre-calibration relative intensities: B=0.796, V=0.204
        # Post-calibration: V = round(0.8401 * 0.204, 3) = 0.171, B = 1 - 0.171 = 0.829
        assert np.isclose(ints[0], 0.0, atol=1e-3)
        assert np.isclose(ints[1], 0.829, atol=1e-3)
        assert np.isclose(ints[2], 0.171, atol=1e-3)
        # Verify pre-calibration ratio recovered from historical CSV
        assert np.isclose(ints[2] / 0.8401, 0.204, atol=2e-3)

    @pytest.mark.skipif(GLOBAL_SOFTWARE_DIR is None or not GLOBAL_SOFTWARE_DIR.exists(), reason="Global software benchmark directory not configured or found")
    def test_p34hb_dconpeaks_golden_match(self, tmp_path):
        """
        Verify that analyzing real P34HB dconpeaks.txt (CF A1000P Mix2 7Nov25)
        produces exact calibrated monomer estimation in PHACS-NMR.
        """
        sample_dir = GLOBAL_SOFTWARE_DIR / "CF A1000P Mix2 7Nov25" / "10" / "pdata" / "1"
        dcon_file = sample_dir / "dconpeaks.txt"
        historical_csv = sample_dir / "CF A1000P Mix2 7Nov25_dconpeaks_Software_results_13-11-2025 08-39-44.csv"

        assert dcon_file.exists(), f"File not found: {dcon_file}"
        assert historical_csv.exists(), f"File not found: {historical_csv}"

        # Copy to isolated tmp_path preserving depth >= 5 components
        nested_tmp = tmp_path / "global" / "sample" / "CF A1000P Mix2 7Nov25" / "10" / "pdata" / "1"
        nested_tmp.mkdir(parents=True)
        test_file = nested_tmp / dcon_file.name
        test_file.write_bytes(dcon_file.read_bytes())

        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(test_file), sensitivity=0.01)

        assert polymer_type == "PH-F-B"
        assert p_id[0] == 1  # 4HB
        assert p_id[1] == 1  # 3HB

        # In 2026 model: 4HB is scaled by 0.5437
        # 4.399 / (4.399 + 2.304) = 0.656 -> 0.5437 * 0.656 = 0.357, B = 0.643
        assert np.isclose(ints[0], 0.357, atol=1e-3)
        assert np.isclose(ints[1], 0.643, atol=1e-3)
