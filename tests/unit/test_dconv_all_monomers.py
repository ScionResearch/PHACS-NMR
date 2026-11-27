from pathlib import Path
import pytest
import numpy as np
from phacs_nmr.core import io as Read_dcon_file


def create_dcon_file(filepath: Path, peak_records: list[tuple[float, float, float]]):
    """
    Helper to write valid TopSpin dconpeaks.txt format.
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


class TestDconvAllMonomers:
    def test_terpolymer_4hb_3hb_3hv(self, tmp_path):
        """Test branch p_ID == [1,1,1,0,0,0,0,0,0] (4HB + 3HB + 3HV)."""
        sample_file = tmp_path / "data" / "run" / "terpolymer1" / "10" / "pdata" / "1" / "dconpeaks.txt"
        peaks = [
            (4.1085, 20000.0, 50000.0),  # F Est
            (2.3520, 20000.0, 50000.0),  # F CH2
            (5.2620, 50000.0, 150000.0), # B Est
            (2.5975, 50000.0, 150000.0), # B CH2
            (5.1530, 30000.0, 90000.0),  # V Est
            (2.5580, 30000.0, 90000.0),  # V CH2
        ]
        create_dcon_file(sample_file, peaks)
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(sample_file), sensitivity=0.01)

        assert p_id[:3] == [1, 1, 1]
        assert sum(p_id[3:]) == 0
        assert polymer_type == "PH-F-B-V"
        # Total fractions must sum to 1.0
        assert np.isclose(ints[0] + ints[1] + ints[2], 1.0, atol=1e-3)

    def test_terpolymer_3hb_3hhx_3ho(self, tmp_path):
        """Test branch p_ID == [0,1,0,1,0,1,0,0,0] (3HB + 3HHx + 3HO)."""
        sample_file = tmp_path / "data" / "run" / "terpolymer2" / "10" / "pdata" / "1" / "dconpeaks.txt"
        peaks = [
            (5.2620, 50000.0, 150000.0), # B Est
            (1.2700, 50000.0, 150000.0),  # B Met
            (5.2050, 25000.0, 75000.0),  # X Est
            (0.9110, 25000.0, 75000.0),   # X Met
            (5.1820, 25000.0, 75000.0),  # O Est
            (0.8790, 25000.0, 75000.0),   # O Met
        ]
        create_dcon_file(sample_file, peaks)
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(sample_file), sensitivity=0.01)

        assert p_id[1] == 1  # B
        assert p_id[3] == 1  # X
        assert p_id[5] == 1  # O
        assert polymer_type == "PH-B-X-O"
        assert np.isclose(ints[1] + ints[3] + ints[5], 1.0, atol=1e-3)

    def test_complex_mcl_monomers_general_branch(self, tmp_path):
        """Test general branch with Hp, N, and >C9 monomers."""
        sample_file = tmp_path / "data" / "run" / "mcl_sample" / "10" / "pdata" / "1" / "dconpeaks.txt"
        peaks = [
            (5.2620, 60000.0, 180000.0), # B Est (provides total_int)
            (1.2700, 60000.0, 180000.0),  # B Met
            (0.9010, 15000.0, 45000.0),   # Hp Met
            (0.8760, 15000.0, 45000.0),   # N Met
            (0.8769, 10000.0, 30000.0),   # >C9 Met
        ]
        create_dcon_file(sample_file, peaks)
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(sample_file), sensitivity=0.01)

        assert p_id[1] == 1  # B
        assert p_id[4] == 1  # Hp
        assert p_id[6] == 1  # N
        assert p_id[7] == 1  # >C9
        assert "PH-B-Hp-N->C9" in polymer_type

    def test_no_monomers_found_fallback(self, tmp_path):
        """Test fallback when peaks do not match any known monomer."""
        sample_file = tmp_path / "data" / "run" / "empty_sample" / "10" / "pdata" / "1" / "dconpeaks.txt"
        peaks = [
            (7.5000, 10000.0, 30000.0),  # Aromatic peak outside all PHA windows
        ]
        create_dcon_file(sample_file, peaks)
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(sample_file), sensitivity=0.01)
        assert p_id == [0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert polymer_type == "PH"
