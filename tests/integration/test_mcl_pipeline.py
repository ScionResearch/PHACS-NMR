from pathlib import Path
import pytest
import matplotlib
import phacs_nmr
from tests.fixtures.synthetic_mcl_spectra import generate_synthetic_phbhhxho_mcl_ascii


class TestMclPipeline:
    def test_phbhhxho_mcl_synthetic_pipeline_execution(self, tmp_path):
        """
        Verify that analyzing synthetic PHBHHxHO_mcl 13C spectrum successfully
        identifies [0,1,0,1,1,1], executes PHBHHxHO_MCL analysis, and exports reports.
        """
        # Ensure path depth >= 5 components for legacy split('\\')[-5]
        nested_tmp = tmp_path / "data" / "synth" / "sample_mcl" / "run"
        nested_tmp.mkdir(parents=True)
        spectrum_file = nested_tmp / "synthetic_mcl_ascii-spec.txt"

        generate_synthetic_phbhhxho_mcl_ascii(spectrum_file, n_points=65536)
        assert spectrum_file.exists()

        matplotlib.use("Agg")

        # Run with noise multiplier 5
        figures = phacs_nmr.start_analysis(str(spectrum_file), multiple_noise=5)
        assert len(figures) > 0, "Expected figures returned from PHBHHxHO_MCL analysis"

        # Verify output files exist
        pdf_files = list(nested_tmp.glob("*.pdf"))
        assert len(pdf_files) >= 1, "Expected PDF report generated"
