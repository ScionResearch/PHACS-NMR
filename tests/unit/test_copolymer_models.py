"""
Unit and invariant tests for copolymer microstructure statistical modeling.
Tests core.models functions directly across PHBV, P34HB, and PHBHHx profiles.
Verifies mathematical invariants (normalization, diad/triad consistency, Bernoullian D=1 and chi=1),
Model 3 blend parameter recovery, and numerical guards against division by zero.
"""
import pytest
import numpy as np

from phacs_nmr.core.models import (
    calculate_d_parameter,
    calculate_chi_parameter,
    calculate_average_lengths,
    compute_kamiya_model_1,
    compute_kamiya_model_2,
    compute_kamiya_model_3,
    apply_copolymer_statistics,
)
from phacs_nmr.core.config import PHBV_PROFILE, P34HB_PROFILE, PHBHHX_PROFILE


class TestCopolymerModels:
    """Test suite for copolymer mathematical models and statistical invariants."""

    # -------------------------------------------------------------------------
    # 1. Parameter calculation unit tests & division-by-zero guards
    # -------------------------------------------------------------------------

    def test_d_parameter_calculation_and_sentinels(self):
        """Verify Kamiya D calculation and sentinel (100.0) when cross-diads are zero or NaN."""
        keys = ("AA", "BB", "AB", "BA")

        # Ideal random: F_AA * F_BB == F_AB * F_BA => D == 1.0
        diads = {"AA": 0.09, "BB": 0.49, "AB": 0.21, "BA": 0.21}
        assert np.isclose(calculate_d_parameter(diads, *keys), 1.0, atol=1e-6)

        # Blocky: higher homodiads => D > 1.0
        blocky = {"AA": 0.25, "BB": 0.65, "AB": 0.05, "BA": 0.05}
        assert calculate_d_parameter(blocky, *keys) > 1.0

        # Alternating: higher cross-diads => D < 1.0
        alternating = {"AA": 0.05, "BB": 0.05, "AB": 0.45, "BA": 0.45}
        assert calculate_d_parameter(alternating, *keys) < 1.0

        # Homopolymer guards: cross-diads are zero => 100.0 sentinel
        homopolymer = {"AA": 1.0, "BB": 0.0, "AB": 0.0, "BA": 0.0}
        assert calculate_d_parameter(homopolymer, *keys) == 100.0

        # Missing or NaN cross-diads => 100.0 sentinel
        nan_diads = {"AA": 0.5, "BB": 0.5, "AB": np.nan, "BA": 0.2}
        assert calculate_d_parameter(nan_diads, *keys) == 100.0

    def test_chi_parameter_calculation_and_sentinels(self):
        """Verify Stanford chi parameter calculation and NaN sentinels."""
        keys = ("A*", "B*", "BA")

        # Ideal random: chi == 1.0
        diads = {"A*": 0.3, "B*": 0.7, "BA": 0.21}
        assert np.isclose(calculate_chi_parameter(diads, *keys), 1.0, atol=1e-6)

        # Alternating: chi > 1.0
        diads_alt = {"A*": 0.5, "B*": 0.5, "BA": 0.45}
        assert calculate_chi_parameter(diads_alt, *keys) > 1.0

        # Blocky: chi < 1.0
        diads_blk = {"A*": 0.5, "B*": 0.5, "BA": 0.05}
        assert calculate_chi_parameter(diads_blk, *keys) < 1.0

        # Zero monomer fraction => returns NaN without division error
        zero_a = {"A*": 0.0, "B*": 1.0, "BA": 0.0}
        assert np.isnan(calculate_chi_parameter(zero_a, *keys))

        # NaN monomer fraction => returns NaN
        nan_diad = {"A*": np.nan, "B*": 0.5, "BA": 0.2}
        assert np.isnan(calculate_chi_parameter(nan_diad, *keys))

    def test_average_lengths_calculation_and_sentinels(self):
        """Verify number-average sequence lengths len_A and len_B."""
        keys = ("A*", "B*", "AB", "BA")

        diads = {"A*": 0.2, "B*": 0.8, "AB": 0.16, "BA": 0.16}
        len_a, len_b = calculate_average_lengths(diads, *keys)
        assert np.isclose(len_a, 0.2 / 0.16, atol=1e-6)
        assert np.isclose(len_b, 0.8 / 0.16, atol=1e-6)

        # Zero cross-diad => returns NaN without crash
        zero_cross = {"A*": 0.2, "B*": 0.8, "AB": 0.0, "BA": 0.0}
        len_a_zero, len_b_zero = calculate_average_lengths(zero_cross, *keys)
        assert np.isnan(len_a_zero) and np.isnan(len_b_zero)

    # -------------------------------------------------------------------------
    # 2. Model 1 (Bernoullian) mathematical invariants across profiles
    # -------------------------------------------------------------------------

    @pytest.mark.parametrize("profile", [PHBV_PROFILE, P34HB_PROFILE, PHBHHX_PROFILE])
    def test_model_1_invariants_across_all_profiles(self, profile):
        """
        Verify that compute_kamiya_model_1 satisfies probability conservation,
        triad-to-diad marginals, and D == 1.0, chi == 1.0 for any composition.
        """
        km = profile.keys_map

        for f_a in [0.05, 0.2, 0.5, 0.85]:
            f_b = 1.0 - f_a
            avg_ratio = {km["A*"]: f_a, km["B*"]: f_b}

            m1 = compute_kamiya_model_1(avg_ratio, km)

            # Monomer conservation
            assert np.isclose(m1[km["A*"]], f_a, atol=1e-10)
            assert np.isclose(m1[km["B*"]], f_b, atol=1e-10)

            # Diad normalization: F_AA + F_AB + F_BA + F_BB == 1.0
            sum_diads = m1[km["AA"]] + m1[km["AB"]] + m1[km["BA"]] + m1[km["BB"]]
            assert np.isclose(sum_diads, 1.0, atol=1e-10)

            # Triad normalization
            sum_triads = (
                m1[km["AAA"]] + m1[km["BAA"]] + m1[km["AAB"]] + m1[km["BAB"]] +
                m1[km["BBB"]] + m1[km["ABB"]] + m1[km["BBA"]] + m1[km["ABA"]]
            )
            assert np.isclose(sum_triads, 1.0, atol=1e-10)

            # Marginal consistency: F_AA = F_AAA + F_AAB
            assert np.isclose(m1[km["AA"]], m1[km["AAA"]] + m1[km["AAB"]], atol=1e-10)
            # F_BA = F_BAA + F_BAB
            assert np.isclose(m1[km["BA"]], m1[km["BAA"]] + m1[km["BAB"]], atol=1e-10)
            # F_BB = F_BBB + F_BBA
            assert np.isclose(m1[km["BB"]], m1[km["BBB"]] + m1[km["BBA"]], atol=1e-10)

            # Invariant: D == 1.0 and chi == 1.0 for Bernoullian random copolymer
            d_val = calculate_d_parameter(m1, km["AA"], km["BB"], km["AB"], km["BA"])
            chi_val = calculate_chi_parameter(m1, km["A*"], km["B*"], km["BA"])
            assert np.isclose(d_val, 1.0, atol=1e-10)
            assert np.isclose(chi_val, 1.0, atol=1e-10)

    # -------------------------------------------------------------------------
    # 3. Model 2 (Markovian 1st Order) invariants & recovery
    # -------------------------------------------------------------------------

    def test_model_2_recovers_model_1_for_random_copolymer(self):
        """When fed random copolymer diads, Model 2 must produce identical results to Model 1."""
        km = PHBV_PROFILE.keys_map
        f_v = 0.35
        f_b = 0.65

        # Feed Bernoullian diads to Model 2
        avg_ratio = {
            km["A*"]: f_v,
            km["B*"]: f_b,
            km["AA"]: f_v ** 2,
            km["BB"]: f_b ** 2,
            km["AB"]: f_v * f_b,
            km["BA"]: f_v * f_b,
        }

        m1 = compute_kamiya_model_1(avg_ratio, km)
        m2 = compute_kamiya_model_2(avg_ratio, km)

        for key in [km["AA"], km["BB"], km["AB"], km["BA"], km["AAA"], km["BBB"], km["A*"], km["B*"]]:
            assert np.isclose(m2[key], m1[key], atol=1e-6), f"Mismatch on {key}: m2={m2[key]}, m1={m1[key]}"

        # Conditional probabilities for ideal random: Paa = f_a, Pbb = f_b
        assert np.isclose(m2[km["Paa"]], f_v, atol=1e-6)
        assert np.isclose(m2[km["Pbb"]], f_b, atol=1e-6)

    def test_model_2_zero_monomer_guard(self):
        """Model 2 with zero monomer fraction should not crash with ZeroDivisionError."""
        km = PHBV_PROFILE.keys_map
        avg_ratio = {
            km["A*"]: 0.0,
            km["B*"]: 1.0,
            km["AA"]: 0.0,
            km["BB"]: 1.0,
            km["AB"]: 0.0,
            km["BA"]: 0.0,
        }
        m2 = compute_kamiya_model_2(avg_ratio, km)
        assert m2[km["Paa"]] == 0.0
        assert m2[km["Pbb"]] == 1.0
        # All output values must be finite numbers without raising ZeroDivisionError
        for k, v in m2.items():
            assert not np.isnan(v), f"Value for {k} is NaN"

    # -------------------------------------------------------------------------
    # 4. Model 3 (Blend) recovery via fsolve
    # -------------------------------------------------------------------------

    def test_model_3_fsolve_synthetic_blend_recovery(self):
        """
        Verify that compute_kamiya_model_3 correctly solves Kamiya Eq (6)
        and recovers known blend parameters A, B, and X.
        """
        km = PHBV_PROFILE.keys_map
        # Known true parameters: Copolymer A (HV=0.7) + Copolymer B (HV=0.2) in ratio X=0.6 : 0.4
        a_true, b_true, x_true = 0.70, 0.20, 0.60

        f_vvv = (a_true**3) * x_true + (b_true**3) * (1.0 - x_true)
        f_bvv = (a_true**2) * (1.0 - a_true) * x_true + (b_true**2) * (1.0 - b_true) * (1.0 - x_true)
        f_bvb = a_true * ((1.0 - a_true)**2) * x_true + b_true * ((1.0 - b_true)**2) * (1.0 - x_true)

        avg_ratio = {
            km["AAA"]: f_vvv,
            km["BAA"]: f_bvv,
            km["AAB"]: f_bvv,
            km["BAB"]: f_bvb,
        }

        m3 = compute_kamiya_model_3(avg_ratio, km)

        # Check solution recovered
        recovered_fractions = sorted([m3["A"], m3["B"]])
        expected_fractions = sorted([a_true, b_true])
        assert np.allclose(recovered_fractions, expected_fractions, atol=1e-3)

        # Check triad reconstruction consistency
        assert np.isclose(m3[km["AAA"]], f_vvv, atol=1e-4)
        assert np.isclose(m3[km["BAB"]], f_bvb, atol=1e-4)

    # -------------------------------------------------------------------------
    # 5. apply_copolymer_statistics integration across table rows
    # -------------------------------------------------------------------------

    def test_apply_copolymer_statistics_rows(self):
        """Verify apply_copolymer_statistics annotates D, chi, len_A, len_B on all rows."""
        km = PHBV_PROFILE.keys_map

        ratio_dic = {
            1: {"V*V": 0.04, "B*B": 0.64, "V*B": 0.16, "B*V": 0.16, "V*": 0.2, "B*": 0.8},
            "avg": {"V*V": 0.04, "B*B": 0.64, "V*B": 0.16, "B*V": 0.16, "V*": 0.2, "B*": 0.8},
            "model_1": {"V*V": 0.04, "B*B": 0.64, "V*B": 0.16, "B*V": 0.16, "V*": 0.2, "B*": 0.8},
        }

        result = apply_copolymer_statistics(ratio_dic, km)

        for row_key in [1, "avg", "model_1"]:
            row = result[row_key]
            assert "D" in row
            assert "chi" in row
            assert "len_V" in row
            assert "len_B" in row
            assert np.isclose(row["D"], 1.0, atol=1e-6)
            assert np.isclose(row["chi"], 1.0, atol=1e-6)
            assert np.isclose(row["len_V"], 0.2 / 0.16, atol=1e-6)
            assert np.isclose(row["len_B"], 0.8 / 0.16, atol=1e-6)
