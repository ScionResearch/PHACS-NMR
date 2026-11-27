"""
Unit tests for MainWindow GUI components, layout, and event handlers.
"""
import sys
import os
from pathlib import Path
import pytest
from PyQt6.QtCore import QUrl, QMimeData, Qt

from phacs_nmr.gui import MainWindow, QtOutputRedirector


class TestMainWindow:
    """Test suite for MainWindow UI behavior, drag-and-drop, and layout management."""

    def test_window_initialization(self, qapp):
        """MainWindow should initialize with title, default widget states, and read-only console."""
        window = MainWindow()
        try:
            assert window.windowTitle() == "PHACS-NMR"
            assert not window.windowIcon().isNull()
            assert window.acceptDrops() is True
            assert window.btn_C_NMR_Analysis.isEnabled() is True
            assert window.btn_H_NMR_Analysis.isEnabled() is True
            assert window.btn_Help.isEnabled() is True
            assert window.progress_bar.isHidden() is True
            assert window.button_back.isHidden() is True
            assert window.pdf_button.isHidden() is True
            assert window.text_edit.isReadOnly() is True
            assert "Ready" in window.status_label.text()
        finally:
            window.close()

    def test_help_toggle_preserves_layout(self, qapp):
        """Toggling Help and Back menu should update visibility without breaking layout hierarchy."""
        window = MainWindow()
        try:
            # Click Help
            window.help_clicked()
            assert window.button_back.isHidden() is False
            assert window.pdf_button.isHidden() is False
            assert "Displaying Help" in window.status_label.text()
            assert len(window.text_edit.toPlainText()) > 0

            # Click Back
            window.back_menu()
            assert window.button_back.isHidden() is True
            assert window.pdf_button.isHidden() is True
            assert "Returned to main log view" in window.text_edit.toPlainText()
            assert window.status_label.text() == "Ready"

            # Click Help again to ensure widgets were not orphaned or deleted from layout
            window.help_clicked()
            assert window.button_back.isHidden() is False
            assert window.pdf_button.isHidden() is False
        finally:
            window.close()

    def test_logging_and_clear_log(self, qapp):
        """Appending to log should update console, and clear_log should empty it."""
        window = MainWindow()
        try:
            window.clear_log()
            assert window.text_edit.toPlainText() == ""

            window.append_log("Test log entry 123")
            assert "Test log entry 123" in window.text_edit.toPlainText()

            window.clear_log()
            assert window.text_edit.toPlainText() == ""
        finally:
            window.close()

    def test_drag_mime_validation(self, qapp, tmp_path):
        """is_valid_drag_mime should accept txt/csv files and reject non-matching files."""
        window = MainWindow()
        try:
            # 1. Valid .txt file
            txt_file = tmp_path / "ascii-spec.txt"
            txt_file.write_text("dummy")
            mime_txt = QMimeData()
            mime_txt.setUrls([QUrl.fromLocalFile(str(txt_file))])
            assert window.is_valid_drag_mime(mime_txt) is True

            # 2. Valid .csv file
            csv_file = tmp_path / "results.csv"
            csv_file.write_text("dummy")
            mime_csv = QMimeData()
            mime_csv.setUrls([QUrl.fromLocalFile(str(csv_file))])
            assert window.is_valid_drag_mime(mime_csv) is True

            # 3. Non-txt file (e.g. .png)
            png_file = tmp_path / "image.png"
            png_file.write_text("dummy")
            mime_png = QMimeData()
            mime_png.setUrls([QUrl.fromLocalFile(str(png_file))])
            assert window.is_valid_drag_mime(mime_png) is False

            # 4. Empty mime data
            mime_empty = QMimeData()
            assert window.is_valid_drag_mime(mime_empty) is False
        finally:
            window.close()

    def test_stdout_and_stderr_redirection_and_cleanup(self, qapp):
        """Writing to sys.stdout or sys.stderr while window is open should pipe to console and restore on close."""
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        window = MainWindow()
        try:
            assert isinstance(sys.stdout, QtOutputRedirector)
            assert isinstance(sys.stderr, QtOutputRedirector)
            print("Intercepted stdout line 456")
            sys.stderr.write("Intercepted stderr warning 789\n")
            qapp.processEvents()
            assert "Intercepted stdout line 456" in window.text_edit.toPlainText()
            assert "Intercepted stderr warning 789" in window.text_edit.toPlainText()
        finally:
            window.close()
            assert sys.stdout is original_stdout
            assert sys.stderr is original_stderr

    def test_phacs_nmr_gui_launcher_import(self):
        """Verify that gui.py can be loaded dynamically and exposes expected attributes."""
        import importlib.util
        launcher_path = Path(__file__).resolve().parent.parent.parent / "src" / "phacs_nmr" / "gui.py"
        assert launcher_path.exists(), f"Launcher not found: {launcher_path}"

        spec = importlib.util.spec_from_file_location("phacs_nmr_gui_module", str(launcher_path))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        assert hasattr(mod, "MainWindow")
        assert hasattr(mod, "AnalysisWorker")
        assert hasattr(mod, "main")

    def test_phacs_nmr_launcher_import(self):
        """Verify that cli.py can be loaded dynamically and exposes expected attributes."""
        import importlib.util
        launcher_path = Path(__file__).resolve().parent.parent.parent / "src" / "phacs_nmr" / "cli.py"
        assert launcher_path.exists(), f"Launcher not found: {launcher_path}"

        spec = importlib.util.spec_from_file_location("phacs_nmr_cli_module", str(launcher_path))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        assert hasattr(mod, "start_analysis")
        assert hasattr(mod, "main")

