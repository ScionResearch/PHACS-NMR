"""
Unit tests for path_utils.py.
Verifies safe extraction of sample names and output prefixes across Bruker TopSpin
directory layouts, non-TopSpin deep and shallow paths, default/custom filenames, and flat files.
"""
from pathlib import Path
from phacs_nmr.path_utils import resolve_sample_name, resolve_output_prefix, is_drive_or_root


class TestPathUtils:
    def test_is_drive_or_root(self):
        assert is_drive_or_root("D:") is True
        assert is_drive_or_root("D:\\") is True
        assert is_drive_or_root("C:/") is True
        assert is_drive_or_root("/") is True
        assert is_drive_or_root("\\") is True
        assert is_drive_or_root("") is True
        assert is_drive_or_root("Sample123") is False
        assert is_drive_or_root("PHBV_software - 13C") is False

    def test_topspin_standard_hierarchy_default_dconpeaks(self):
        """TopSpin: <root>/<sample_name>/<expno>/pdata/<procno>/dconpeaks.txt"""
        p = Path(r"C:\Users\FornesC\NMR\Commercial samples\CF BP 350-PD PHBHHx Aug25\10\pdata\1\dconpeaks.txt")
        assert resolve_sample_name(p) == "CF BP 350-PD PHBHHx Aug25"
        # Default filename 'dconpeaks' is omitted from output prefix
        assert resolve_output_prefix(p) == "CF BP 350-PD PHBHHx Aug25"

    def test_topspin_forward_slash_default_ascii_spec(self):
        """POSIX or forward slash paths with pdata hierarchy and default ascii-spec.txt."""
        p = "C:/Users/FornesC/NMR/Commercial samples/250911_CF_UQ_PHBV8/11/pdata/1/ascii-spec.txt"
        assert resolve_sample_name(p) == "250911_CF_UQ_PHBV8"
        # Default filename 'ascii-spec' is omitted from output prefix
        assert resolve_output_prefix(p) == "250911_CF_UQ_PHBV8"

    def test_topspin_custom_filestem_included(self):
        """TopSpin hierarchy with custom non-default file stem must include the stem."""
        p = Path(r"C:\Users\FornesC\NMR\Commercial samples\CF BP 350-PD PHBHHx Aug25\10\pdata\1\carbon13_proc.txt")
        assert resolve_sample_name(p) == "CF BP 350-PD PHBHHx Aug25"
        assert resolve_output_prefix(p) == "CF BP 350-PD PHBHHx Aug25_carbon13_proc"

    def test_topspin_dconpeak_singular_treated_as_custom(self):
        """dconpeak.txt (singular) is NOT a default stem; it must be included in output prefix."""
        p = Path(r"C:\Users\FornesC\NMR\Commercial samples\CF BP 350-PD PHBHHx Aug25\10\pdata\1\dconpeak.txt")
        assert resolve_sample_name(p) == "CF BP 350-PD PHBHHx Aug25"
        assert resolve_output_prefix(p) == "CF BP 350-PD PHBHHx Aug25_dconpeak"

    def test_no_pdata_deep_path_uses_filestem(self):
        """When pdata is absent in a deep path (depth >= 5), must use filestem, NOT parts[-5]."""
        p = Path(r"C:\Users\jsmith\Documents\Projects\PHA_Research\my_sample_run.txt")
        assert resolve_sample_name(p) == "my_sample_run"
        assert resolve_output_prefix(p) == "my_sample_run"

    def test_flat_path_avoid_drive_letter(self):
        """Path with 5 parts where parts[-5] is a drive letter (e.g. D:\\mock_data\\benchmarks\\PHBV\\sample.txt)."""
        p = Path(r"D:\mock_data\benchmarks\PHBV\sample.txt")
        name = resolve_sample_name(p)
        assert not is_drive_or_root(name), f"Returned drive letter or root: {name}"
        assert name == "sample"
        assert resolve_output_prefix(p) == "sample"

    def test_shallow_path_no_index_error(self):
        """Path with depth < 5 must never raise IndexError and returns filestem when no pdata."""
        p = Path(r"C:\temp\sample.txt")
        name = resolve_sample_name(p)
        assert not is_drive_or_root(name)
        assert name == "sample"
        assert resolve_output_prefix(p) == "sample"

    def test_single_file_fallback(self):
        """Bare file in current directory with no pdata returns filestem."""
        p = Path("ascii-spec.txt")
        name = resolve_sample_name(p)
        assert name == "ascii-spec"
        assert resolve_output_prefix(p) == "ascii-spec"
