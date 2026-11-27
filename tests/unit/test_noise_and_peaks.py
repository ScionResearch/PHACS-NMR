import pytest
import numpy as np
from pathlib import Path
import phacs_nmr as dnd


class TestNoiseAndPeaks:
    def test_inv_scale_exact_and_closest(self):
        """inv_scale should find closest index in array."""
        arr = np.array([10.0, 8.0, 6.0, 4.0, 2.0, 0.0])
        assert dnd.inv_scale(arr, 8.0) == 1
        assert dnd.inv_scale(arr, 7.9) == 1
        assert dnd.inv_scale(arr, 6.1) == 2
        assert dnd.inv_scale(arr, 15.0) == 0
        assert dnd.inv_scale(arr, -5.0) == 5

    def test_calc_noise_level(self):
        """calc_noise_level should estimate noise in region -5 to -15 ppm."""
        ppm_scale = np.linspace(200.0, -20.0, 2200)
        # Create pure zero data with known noise in -5 to -15 region
        data = np.zeros_like(ppm_scale)
        mask = (ppm_scale <= -5.0) & (ppm_scale >= -15.0)
        np.random.seed(42)
        noise_samples = np.random.normal(0, 1000.0, mask.sum())
        data[mask] = noise_samples

        noise = dnd.calc_noise_level(ppm_scale, data)
        # Expected mean absolute value for normal distribution is sigma * sqrt(2/pi) ~ 797.88
        assert 700.0 < noise < 900.0

    def test_treatment_of_txt_file(self, tmp_path):
        """treatment_of_txt_file should parse 4-column ASCII TopSpin spectrum."""
        test_file = tmp_path / "ascii-spec.txt"
        lines = [
            "1, 1234.5, 0.0, 200.0\n",
            "2, 2345.6, 0.0, 199.5\n",
            "3, 3456.7, 0.0, 199.0\n",
        ]
        test_file.write_text("".join(lines), encoding="cp1252")

        x, y = dnd.treatment_of_txt_file(str(test_file))
        assert len(x) == 3
        assert len(y) == 3
        # PPM in col 3 (index 3)
        assert x == [200.0, 199.5, 199.0]
        # Intensity in col 1 (index 1)
        assert y == [1234.5, 2345.6, 3456.7]
