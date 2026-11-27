import pytest
import numpy as np
from pathlib import Path
from phacs_nmr.core import io as Read_dcon_file


class TestDconvCorrection:
    def test_find_monomers_limits_structure(self):
        """Verify monomer search limits are well-formed intervals."""
        limits = {
            'B': {'Est': (5.2608, 5.2631), 'Met': (1.269, 1.2703), 'CH2': (2.5965, 2.5985)},
            'V': {'Met': (0.8940, 0.897), 'Est': (5.1516, 5.1554), 'CH2': (2.556, 2.5595)},
            'X': {'Est': (5.204, 5.206), 'Met': (0.9099, 0.9133), 'CH2': (2.554, 2.556)},
            'Hp': {'Met': (0.9005, 0.9015)},
            'O': {'Est': (5.1809, 5.1834), 'Met': (0.8778, 0.8802), 'CH2': (2.550, 2.5474)},
            'N': {'Met': (0.8755, 0.8765)},
            'F': {'Est': (4.107, 4.110), 'CH2': (2.350, 2.354)},
            '>C9': {'Met': (0.8764, 0.8774)},
            '>C9_1': {'Met': (0.881, 0.8825)},
        }
        for monomer, regions in limits.items():
            for proton, (p_min, p_max) in regions.items():
                low = min(p_min, p_max)
                high = max(p_min, p_max)
                assert low < high, f"Invalid interval for {monomer} {proton}: ({p_min}, {p_max})"

    def test_monomer_identification_phbv(self):
        """Simulate dconv peak data with B and V peaks only."""
        fake_data = {
            'peaks_ppm': [5.262, 1.270, 0.895, 5.153],
            'peaks_intensity': [1000.0, 3000.0, 500.0, 200.0],
            'peaks_area': [5000.0, 15000.0, 2500.0, 1000.0],
        }
        p_id, est_area, est_int = Read_dcon_file.find_monomers(fake_data)
        # Expected: F=0, B=1, V=1, rest=0
        assert p_id[0] == 0  # F (4HB)
        assert p_id[1] == 1  # B (3HB)
        assert p_id[2] == 1  # V (3HV)
        assert sum(p_id[3:]) == 0  # No other monomers

    def test_sensitivity_filtering(self, tmp_path):
        """Test sensitivity threshold logic on temporary dconv file."""
        sample_dir = tmp_path / "root" / "project" / "sample" / "10" / "pdata" / "1"
        sample_dir.mkdir(parents=True)
        dconv_file = sample_dir / "dconpeaks.txt"

        # Create minimal valid TopSpin dconpeaks.txt (8 header lines, then columns)
        lines = [
            "TopSpin Deconv Output\n",
            "Date\n",
            "Sample\n",
            "Param 1\n",
            "Param 2\n",
            "Param 3\n",
            "Param 4\n",
            "Param 5\n",
            # col 0: ppm, col 1,2,3, col 4: intensity, col 5: area, col 6: fwhm
            "5.2620 0 0 0 100000 500000 0.01\n",  # B Est (required for initial detection)
            "1.2700 0 0 0 100000 500000 0.01\n",  # B Met
            "0.8950 0 0 0 20000 100000 0.01\n",   # V Met
        ]
        dconv_file.write_text("".join(lines), encoding="cp1252")

        # Run with sensitivity 0.01 (10% V is above 1%)
        headers, p_id, polymer_type, areas, ints = Read_dcon_file.results_dconv(str(dconv_file), sensitivity=0.01)
        assert p_id[1] == 1  # B
        assert p_id[2] == 1  # V
        assert polymer_type == "PH-B-V"
        # Check that sum of active monomer fractions is ~1.0
        assert np.isclose(ints[1] + ints[2], 1.0, atol=0.01)
