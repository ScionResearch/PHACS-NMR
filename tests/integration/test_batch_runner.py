import sys
from pathlib import Path
import pytest
import matplotlib
import phacs_nmr
from tests.conftest import BENCHMARK_DATA_DIR as BENCHMARK_DIR


class TestBatchRunner:
    def test_batch_runner_invalid_dir(self, monkeypatch, tmp_path):
        """Verify error is raised when folder does not exist."""
        non_existent = tmp_path / "does_not_exist"
        monkeypatch.setattr(sys, "argv", ["phacs_nmr.py", str(non_existent), "10"])
        with pytest.raises(FileNotFoundError):
            phacs_nmr.main()

    def test_batch_runner_empty_dir(self, monkeypatch, tmp_path, capsys):
        """Verify clean exit when no .txt files are found in directory."""
        empty_dir = tmp_path / "empty_folder"
        empty_dir.mkdir()
        monkeypatch.setattr(sys, "argv", ["phacs_nmr.py", str(empty_dir), "10"])
        phacs_nmr.main()
        captured = capsys.readouterr()
        assert "No .txt files found" in captured.out

    @pytest.mark.skipif(BENCHMARK_DIR is None or not BENCHMARK_DIR.exists(), reason="Benchmark dataset directory not configured or found")
    def test_batch_runner_fast_subset_with_error_isolation(self, monkeypatch, tmp_path, capsys):
        """
        Verify batch runner executes over a representative fast subset (PHBV + Homopolymer + Corrupt file)
        and isolates errors without terminating the batch.
        """
        batch_dir = tmp_path / "data" / "batch" / "run_folder" / "10" / "pdata"
        batch_dir.mkdir(parents=True)

        # 1. Valid PHBV spectrum
        phbv_src = BENCHMARK_DIR / "SI_13Cq_PHA HPH18 25 PHBV_11.txt"
        assert phbv_src.exists()
        (batch_dir / "sample1_phbv.txt").write_bytes(phbv_src.read_bytes())

        # 2. Valid Homopolymer (clean early exit)
        homo_src = BENCHMARK_DIR / "SI_13Cq_PHA Y1000P PHBV_11.txt"
        assert homo_src.exists()
        (batch_dir / "sample2_homo.txt").write_bytes(homo_src.read_bytes())

        # 3. Corrupt spectrum to verify error handling
        (batch_dir / "sample3_corrupt.txt").write_text("CORRUPT DATA INVALID\n", encoding="cp1252")

        matplotlib.use("Agg")
        monkeypatch.setattr(sys, "argv", ["phacs_nmr.py", str(batch_dir), "10"])

        # Run batch analysis
        phacs_nmr.main()

        captured = capsys.readouterr()
        assert "=== Processing sample1_phbv.txt ===" in captured.out
        assert "=== Processing sample2_homo.txt ===" in captured.out
        assert "=== Processing sample3_corrupt.txt ===" in captured.out
        assert "Failed on" in captured.out

        # Verify output files for sample1 were created
        generated_csvs = list(batch_dir.glob("*.csv"))
        assert len(generated_csvs) >= 1
