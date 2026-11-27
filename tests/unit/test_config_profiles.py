"""
Unit tests for core.config declarative PolymerProfile configurations and limit builders.
"""
import pytest
import numpy as np

from phacs_nmr.core.config import (
    get_profile,
    PROFILES,
    PHBV_PROFILE,
    P34HB_PROFILE,
    PHBHHX_PROFILE,
    PHBHHXHO_MCL_PROFILE,
    limits_phbv,
    limits_p34hb,
    limits_phbhhx,
    limits_phbhhxho_mcl,
)


class TestConfigProfiles:
    """Test suite for declarative polymer profiles and limit builders."""

    def test_get_profile_lookup(self):
        """get_profile should resolve case-insensitive polymer profile names."""
        assert get_profile("PHBV") == PHBV_PROFILE
        assert get_profile("phbv") == PHBV_PROFILE
        assert get_profile("P34HB") == P34HB_PROFILE
        assert get_profile("p34hb") == P34HB_PROFILE
        assert get_profile("PHBHHx") == PHBHHX_PROFILE
        assert get_profile("phbhhx") == PHBHHX_PROFILE
        assert get_profile("PHBHHXHO_MCL") == PHBHHXHO_MCL_PROFILE
        assert get_profile("phbhhxho_mcl") == PHBHHXHO_MCL_PROFILE

        with pytest.raises(KeyError, match="Unknown polymer profile"):
            get_profile("UNKNOWN_POLYMER_TYPE")

    def test_all_limit_builders(self):
        """Verify that all limit builders execute and return valid limits_s and limits."""
        ppm_scale = np.linspace(200.0, 0.0, 10000)

        for builder in [limits_phbv, limits_p34hb, limits_phbhhx, limits_phbhhxho_mcl]:
            limits_s, limits = builder(ppm_scale)
            assert len(limits_s) > 0
            assert len(limits) > 0
            # Each entry in limits_s has [title, name, start_idx, end_idx]
            for entry in limits_s:
                assert len(entry) == 4
                assert entry[2] >= 0
                assert entry[3] >= 0
            # Each entry in limits is (end_idx, start_idx)
            for region_name, (idx1, idx2) in limits.items():
                assert idx1 >= 0
                assert idx2 >= 0
