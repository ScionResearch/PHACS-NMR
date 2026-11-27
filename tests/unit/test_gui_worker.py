"""
Unit and concurrency tests for AnalysisWorker (QThread).
Verifies background execution, signal emission, and error isolation in headless mode.
"""
import os
import time
from pathlib import Path
import pytest
from PyQt6.QtCore import QEventLoop, QTimer

from phacs_nmr.gui import AnalysisWorker


class TestAnalysisWorker:
    """Test suite for AnalysisWorker thread lifecycle and signals."""

    def test_worker_file_not_found(self, qapp):
        """Worker should emit failed signal when input file does not exist."""
        worker = AnalysisWorker(
            task_type="13C",
            file_path="non_existent_file_path_12345.txt",
            params={"multiple_noise": 10}
        )

        progress_msgs = []
        finished_results = []
        error_msgs = []

        worker.progress.connect(progress_msgs.append)
        worker.finished.connect(finished_results.append)
        worker.failed.connect(error_msgs.append)

        loop = QEventLoop()
        worker.finished.connect(lambda _: loop.quit())
        worker.failed.connect(lambda _: loop.quit())

        worker.start()
        loop.exec()
        worker.wait()

        assert len(error_msgs) == 1
        assert "does not exist" in error_msgs[0]
        assert len(finished_results) == 0

    def test_worker_invalid_task_type(self, qapp, tmp_path):
        """Worker should emit failed signal when task_type is unrecognized."""
        dummy_file = tmp_path / "dummy.txt"
        dummy_file.write_text("dummy content")

        worker = AnalysisWorker(
            task_type="UNKNOWN_TYPE",
            file_path=str(dummy_file),
            params={}
        )

        error_msgs = []
        worker.failed.connect(error_msgs.append)

        loop = QEventLoop()
        worker.finished.connect(lambda _: loop.quit())
        worker.failed.connect(lambda _: loop.quit())

        worker.start()
        loop.exec()
        worker.wait()

        assert len(error_msgs) == 1
        assert "Unknown task type" in error_msgs[0]

    def test_worker_1h_deconvolution_success(self, qapp, tmp_path):
        """Worker should execute 1H deconvolution and emit finished signal."""
        # Create minimal synthetic dconpeaks file
        dcon_file = tmp_path / "sample_dconpeaks.txt"
        dcon_content = """# Bruker deconv peaks
1.278   5.000   100.0   0.05
4.150   10.000  200.0   0.05
"""
        dcon_file.write_text(dcon_content)

        worker = AnalysisWorker(
            task_type="1H",
            file_path=str(dcon_file),
            params={"sensitivity": 0.01}
        )

        progress_msgs = []
        finished_results = []
        error_msgs = []

        worker.progress.connect(progress_msgs.append)
        worker.finished.connect(finished_results.append)
        worker.failed.connect(error_msgs.append)

        loop = QEventLoop()
        worker.finished.connect(lambda _: loop.quit())
        worker.failed.connect(lambda _: loop.quit())

        worker.start()
        loop.exec()
        worker.wait()

        assert len(error_msgs) == 0
        assert len(finished_results) == 1
        res = finished_results[0]
        assert res["task_type"] == "1H"
        assert res["file_path"] == str(dcon_file)
        assert any("1H NMR deconvolution finished successfully" in msg for msg in progress_msgs)

    def test_worker_13c_analysis_success(self, qapp, benchmark_data_dir, tmp_path):
        """Worker should execute 13C analysis in background thread and emit finished signal."""
        # Use verified real benchmark spectrum
        input_txt = benchmark_data_dir / "SI_13Cq_PHA HPH18 25 PHBV_11.txt"
        if not input_txt.exists():
            pytest.skip(f"Benchmark file not found: {input_txt}")

        test_file = tmp_path / input_txt.name
        test_file.write_bytes(input_txt.read_bytes())

        worker = AnalysisWorker(
            task_type="13C",
            file_path=str(test_file),
            params={"multiple_noise": 10}
        )

        progress_msgs = []
        finished_results = []
        error_msgs = []

        worker.progress.connect(progress_msgs.append)
        worker.finished.connect(finished_results.append)
        worker.failed.connect(error_msgs.append)

        loop = QEventLoop()
        worker.finished.connect(lambda _: loop.quit())
        worker.failed.connect(lambda _: loop.quit())

        worker.start()
        loop.exec()
        worker.wait()

        assert len(error_msgs) == 0
        assert len(finished_results) == 1
        res = finished_results[0]
        assert res["task_type"] == "13C"
        assert res["file_path"] == str(test_file)
        assert res["figures_count"] >= 4
        assert any("13C NMR analysis finished successfully" in msg for msg in progress_msgs)
