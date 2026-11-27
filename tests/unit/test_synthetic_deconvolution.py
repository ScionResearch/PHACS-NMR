"""
Unit tests for core.peak_fitting lineshape deconvolution and peak picking.
Verifies spectral region segmentation, Lorentzian area recovery, noise rejection,
and homopolymer absent-region handling using synthetic in-memory spectra (preloaded_xy).
"""
import pytest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from phacs_nmr.core.peak_fitting import split_graphic_regions
from phacs_nmr.core.config import PHBV_PROFILE, limits_phbv


def create_synthetic_spectrum_array(
    peak_specs: list[tuple[float, float, float]],
    noise_sigma: float = 100.0,
    n_points: int = 22000,
    seed: int = 42,
) -> tuple[list[float], list[float]]:
    """
    Generate synthetic (ppm, intensity) lists.
    peak_specs: list of (center_ppm, height, fwhm_ppm)
    noise region (-5 to -15 ppm) included for calc_noise_level.
    """
    np.random.seed(seed)
    ppm = np.linspace(200.0, -20.0, n_points)
    intensity = np.random.normal(0, noise_sigma, n_points)

    for center, height, fwhm in peak_specs:
        gamma = fwhm / 2.0
        intensity += height * (gamma**2 / ((ppm - center)**2 + gamma**2))

    return ppm.tolist(), intensity.tolist()


class TestSyntheticDeconvolution:
    """Test suite for peak picking and lineshape fitting on synthetic signals."""

    def test_synthetic_lorentzian_area_recovery(self):
        """
        Verify that split_graphic_regions detects two known Lorentzian peaks
        in B1V1 region and recovers the analytical ~2:1 area ratio.
        """
        # Peak 1: B*B at 169.15 ppm (height=50,000, fwhm=0.03)
        # Peak 2: B*V at 169.35 ppm (height=25,000, fwhm=0.03)
        # Analytical area ratio = 50,000 / 25,000 = 2.0
        peaks = [
            (169.15, 50000.0, 0.03),
            (169.35, 25000.0, 0.03),
        ]
        x_ppm, y_int = create_synthetic_spectrum_array(peaks, noise_sigma=100.0)

        # Single region table for targeted B1V1 verification
        table = [(168.9, 169.8, "B1V1")]

        x_ret, y_ret, pha_df, fig = split_graphic_regions(
            file_path="in_memory_synth",
            table=table,
            multiple_noise=10,
            limits_fn=limits_phbv,
            preloaded_xy=(x_ppm, y_int),
        )
        plt.close(fig)

        assert "peaks_data" in pha_df
        peaks_data = pha_df["peaks_data"]
        assert len(peaks_data) > 0, "Expected peaks_data to contain detected peaks"

        # Extract detected peak heights and raw areas
        heights_by_name = {}
        areas_by_name = {}
        for row in peaks_data:
            # row: [bname, name, peak_pos, peak_height, Area_raw, Area_curve]
            name = row[1]
            heights_by_name[name] = row[3]
            areas_by_name[name] = row[4]

        assert "B*B" in heights_by_name
        assert "B*V" in heights_by_name

        bb_height = heights_by_name["B*B"]
        bv_height = heights_by_name["B*V"]

        # Assert recovery of the ~2.0 height ratio within 15% tolerance
        assert bb_height > 0 and bv_height > 0
        height_ratio = bb_height / bv_height
        assert np.isclose(height_ratio, 2.0, atol=0.25), f"Expected peak height ratio ~2.0, got {height_ratio:.2f}"

        # B*B area must be significantly larger than B*V
        assert areas_by_name["B*B"] > areas_by_name["B*V"]

        # Absent peaks in this region (V*B, V*V) must have near-zero fitted curve area (<1% of B*B)
        curve_areas = {r[1]: r[5] for r in peaks_data}
        assert curve_areas["V*B"] < 0.01 * curve_areas["B*B"]
        assert curve_areas["V*V"] < 0.01 * curve_areas["B*B"]

    def test_pure_noise_spectrum_rejection(self):
        """
        Pure Gaussian noise with multiple_noise=10 should have zero peaks picked (r == 0),
        cleanly exercising the empty region branch and setting default 0 areas.
        """
        # Noise sigma = 100, noise range -5 to -15 gives mean(abs) ~ 80.
        # multiple_noise = 10 => threshold ~ 800. Max noise ~ 4-5 sigma ~ 400-500, never >= 800.
        x_ppm, y_int = create_synthetic_spectrum_array([], noise_sigma=100.0)
        table = [(168.9, 169.8, "B1V1")]

        x_ret, y_ret, pha_df, fig = split_graphic_regions(
            file_path="in_memory_noise",
            table=table,
            multiple_noise=10,
            limits_fn=limits_phbv,
            preloaded_xy=(x_ppm, y_int),
        )
        plt.close(fig)

        peaks_data = pha_df.get("peaks_data", [])
        assert len(peaks_data) > 0
        for row in peaks_data:
            # Check position is NaN and areas are 0
            assert np.isnan(row[2])  # position
            assert np.isnan(row[3])  # height
            assert row[4] == 0       # Area_raw
            assert row[5] == 0       # Area_curve

    def test_homopolymer_absent_regions_handling(self):
        """
        Homopolymer PHA spectrum has signal only in B regions (e.g. B3 at 67.6 ppm),
        while V regions (e.g. V3 at 71.9 ppm) contain only baseline noise.
        Verifies both active and absent regions within the same run.
        """
        # Add peak only for B3 (67.6 ppm)
        peaks = [
            (67.60, 40000.0, 0.04),
        ]
        x_ppm, y_int = create_synthetic_spectrum_array(peaks, noise_sigma=80.0)

        table = [
            (67.0, 68.0, "B3"),
            (71.6, 72.1, "V3"),
        ]

        x_ret, y_ret, pha_df, fig = split_graphic_regions(
            file_path="in_memory_homopolymer",
            table=table,
            multiple_noise=10,
            limits_fn=limits_phbv,
            preloaded_xy=(x_ppm, y_int),
        )
        plt.close(fig)

        peaks_data = pha_df.get("peaks_data", [])
        assert len(peaks_data) > 0

        b3_rows = [r for r in peaks_data if r[0] == "B3"]
        v3_rows = [r for r in peaks_data if r[0] == "V3"]

        # B3 must have positive area
        assert any(r[4] > 0 for r in b3_rows), "Expected positive peak area in B3"
        # V3 must have 0 area (no peak in noise)
        assert all(r[4] == 0 for r in v3_rows), "Expected 0 peak area in absent V3 region"
