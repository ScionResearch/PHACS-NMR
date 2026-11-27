"""
Modernized GUI and concurrency engine for PHACS-NMR.
Provides responsive, thread-safe spectral analysis execution via QThread,
native file dialogs, Drag & Drop support, and integrated reporting
using native platform styling.

Copyright (C) 2026 Bioeconomy Science Institute
Licensed under the GNU General Public License v3.0 or later (GPLv3+).
See LICENSE for details.
"""
import logging
import os
import sys
from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
from PyQt6.QtCore import QObject, Qt, QThread, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices, QDragEnterEvent, QDropEvent, QFont, QIcon
from PyQt6.QtWidgets import (
    QApplication,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from phacs_nmr import cli as phacs_nmr
from phacs_nmr.core import io as dcon_io

logger = logging.getLogger("phacs_nmr.gui")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))


def _find_resource(filename: str) -> str:
    """Find a resource file in package assets or repository root."""
    in_pkg_assets = os.path.join(SCRIPT_DIR, "assets", filename)
    if os.path.isfile(in_pkg_assets):
        return in_pkg_assets
    in_pkg = os.path.join(SCRIPT_DIR, filename)
    if os.path.isfile(in_pkg):
        return in_pkg
    in_root_assets = os.path.join(PROJECT_ROOT, "assets", filename)
    if os.path.isfile(in_root_assets):
        return in_root_assets
    in_root = os.path.join(PROJECT_ROOT, filename)
    if os.path.isfile(in_root):
        return in_root
    return in_pkg_assets


# =============================================================================
# Thread-Safe Output & Logging Interceptor
# =============================================================================

class QtOutputRedirector(QObject):
    """
    Thread-safe stream interceptor redirecting stdout/stderr or logging
    to a Qt signal with line buffering.
    Option B: By default, does not echo to terminal stream, keeping the launching
    terminal clean while the GUI is active.
    """
    text_written = pyqtSignal(str)

    def __init__(self, original_stream=None, echo_to_stream: bool = False):
        super().__init__()
        self.original_stream = original_stream
        self.echo_to_stream = echo_to_stream
        self._buffer = ""

    def write(self, text: str):
        if not text:
            return
        if self.echo_to_stream and self.original_stream:
            try:
                self.original_stream.write(text)
            except Exception:
                pass
        self._buffer += text
        if "\n" in self._buffer:
            lines = self._buffer.split("\n")
            for line in lines[:-1]:
                self.text_written.emit(line)
            self._buffer = lines[-1]

    def flush(self):
        if self._buffer:
            self.text_written.emit(self._buffer)
            self._buffer = ""
        if self.echo_to_stream and self.original_stream:
            try:
                self.original_stream.flush()
            except Exception:
                pass



# =============================================================================
# Background Worker Thread (QThread)
# =============================================================================

class AnalysisWorker(QThread):
    """
    Background worker executing heavy spectral fitting and deconvolution
    without blocking the PyQt main thread.
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    failed = pyqtSignal(str)

    def __init__(self, task_type: str, file_path: str, params: dict[str, Any] | None = None, parent: QObject | None = None):
        super().__init__(parent)
        self.task_type = task_type  # '13C' or '1H'
        self.file_path = str(file_path)
        self.params = params or {}
        self._is_cancelled = False

    def run(self):
        try:
            import matplotlib
            matplotlib.use("Agg")

            if not os.path.isfile(self.file_path):
                raise FileNotFoundError(f"Input file does not exist: {self.file_path}")

            if self.task_type == "13C":
                noise_multiplier = int(self.params.get("multiple_noise", 10))
                self.progress.emit(f"Starting 13C NMR analysis for {os.path.basename(self.file_path)}...")
                self.progress.emit(f"Noise multiplier: {noise_multiplier}")
                
                # Execute analysis
                figures = phacs_nmr.start_analysis(self.file_path, noise_multiplier)
                
                # Find output PDF if generated
                file_dir = os.path.dirname(self.file_path)
                pdf_files = sorted(
                    Path(file_dir).glob("*_Software_results_*.pdf"),
                    key=os.path.getmtime,
                    reverse=True
                )
                pdf_path = str(pdf_files[0]) if pdf_files else None

                # Clean up matplotlib figures to prevent memory leaks and UI collisions
                plt.close("all")

                self.progress.emit("13C NMR analysis finished successfully.")
                self.finished.emit({
                    "task_type": "13C",
                    "file_path": self.file_path,
                    "pdf_path": pdf_path,
                    "figures_count": len(figures) if figures else 0
                })

            elif self.task_type == "1H":
                sensitivity = float(self.params.get("sensitivity", 0.01))
                self.progress.emit(f"Starting 1H NMR deconvolution for {os.path.basename(self.file_path)}...")
                self.progress.emit(f"Sensitivity threshold: {sensitivity}")

                headers, p_id, polymer_type, est_monomers_area, est_monomers_int = (
                    dcon_io.results_dconv(self.file_path, sensitivity)
                )
                self.progress.emit("1H NMR deconvolution finished successfully.")
                self.finished.emit({
                    "task_type": "1H",
                    "file_path": self.file_path,
                    "polymer_type": polymer_type,
                    "p_id": p_id,
                    "est_monomers_int": est_monomers_int,
                })
            else:
                raise ValueError(f"Unknown task type: {self.task_type}")

        except Exception as exc:
            self.progress.emit(f"Error during analysis: {exc}")
            self.failed.emit(str(exc))


# =============================================================================
# Main Window
# =============================================================================

class MainWindow(QWidget):
    """
    Modernized PHACS-NMR application window with non-blocking concurrency,
    drag-and-drop, native controls, and integrated reporting.
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PHACS-NMR")
        self.setMinimumSize(950, 560)
        self.resize(1100, 580)
        self.setAcceptDrops(True)

        # Set application icon
        icon_path = _find_resource("phacs_nmr.ico")
        if os.path.isfile(icon_path):
            self.setWindowIcon(QIcon(icon_path))

        self.worker: AnalysisWorker | None = None
        self.latest_pdf_path: str | None = None
        self.last_directory: str = os.getcwd()

        self.pdf_protocol_path = _find_resource("Software Protocol v2.pdf")
        self.readme_path = _find_resource("READ ME.txt")

        # Set up stdout and stderr redirection (Option B: route exclusively to GUI window, keeping terminal quiet)
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

        self.stdout_redirector = QtOutputRedirector(original_stream=self.original_stdout, echo_to_stream=False)
        self.stdout_redirector.text_written.connect(self.append_log)

        self.stderr_redirector = QtOutputRedirector(original_stream=self.original_stderr, echo_to_stream=False)
        self.stderr_redirector.text_written.connect(self.append_log)

        sys.stdout = self.stdout_redirector
        sys.stderr = self.stderr_redirector

        self._build_ui()

    def _build_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(15, 15, 15, 15)
        main_layout.setSpacing(15)

        # ---------------------------------------------------------
        # Left Panel: Controls and Action Buttons
        # ---------------------------------------------------------
        left_panel = QVBoxLayout()
        left_panel.setSpacing(10)

        # App Brand Title
        title_label = QLabel("PHACS-NMR")
        title_font = QFont("Segoe UI", 16, QFont.Weight.Bold)
        title_label.setFont(title_font)
        left_panel.addWidget(title_label)

        desc_label = QLabel("PHA Polymer Characterization\nby NMR")
        desc_font = QFont("Segoe UI", 9)
        desc_label.setFont(desc_font)
        left_panel.addWidget(desc_label)

        # Separator line
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        left_panel.addWidget(line)

        # Primary Action Buttons
        self.btn_C_NMR_Analysis = QPushButton("13C NMR Analysis")
        self.btn_C_NMR_Analysis.setToolTip("Process 13C ASCII spectrum (ascii-spec.txt) for sequence and microstructure modeling")
        btn_font = QFont()
        btn_font.setBold(True)
        self.btn_C_NMR_Analysis.setFont(btn_font)
        self.btn_C_NMR_Analysis.clicked.connect(self.C_NMR_Analysis_click)
        left_panel.addWidget(self.btn_C_NMR_Analysis)

        self.btn_H_NMR_Analysis = QPushButton("1H NMR Analysis")
        self.btn_H_NMR_Analysis.setToolTip("Process 1H deconvolution peak file (dconpeaks.txt) for monomer composition")
        self.btn_H_NMR_Analysis.setFont(btn_font)
        self.btn_H_NMR_Analysis.clicked.connect(self.H_NMR_Analysis_click)
        left_panel.addWidget(self.btn_H_NMR_Analysis)

        # Open PDF Result Button (initially hidden until analysis finishes)
        self.btn_open_result_pdf = QPushButton("View Generated Report (.pdf)")
        self.btn_open_result_pdf.setFont(btn_font)
        self.btn_open_result_pdf.setVisible(False)
        self.btn_open_result_pdf.clicked.connect(self.open_result_pdf)
        left_panel.addWidget(self.btn_open_result_pdf)

        left_panel.addSpacing(10)

        # Drag and Drop info box
        self.dnd_hint = QLabel("Drag & drop spectrum file\n(.txt / .csv) anywhere on window")
        self.dnd_hint.setFrameShape(QFrame.Shape.Box)
        self.dnd_hint.setFrameShadow(QFrame.Shadow.Sunken)
        self.dnd_hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.dnd_hint.setMargin(10)
        left_panel.addWidget(self.dnd_hint)

        left_panel.addStretch()

        # Help and Protocol Controls
        self.btn_Help = QPushButton("Help & Protocol")
        self.btn_Help.clicked.connect(self.help_clicked)
        left_panel.addWidget(self.btn_Help)

        self.pdf_button = QPushButton("Open Protocol PDF")
        self.pdf_button.clicked.connect(self.open_pdf)
        self.pdf_button.setVisible(False)
        left_panel.addWidget(self.pdf_button)

        self.button_back = QPushButton("<< Back to Log")
        self.button_back.clicked.connect(self.back_menu)
        self.button_back.setVisible(False)
        left_panel.addWidget(self.button_back)

        self.btn_clear_log = QPushButton("Clear Log")
        self.btn_clear_log.clicked.connect(self.clear_log)
        left_panel.addWidget(self.btn_clear_log)

        main_layout.addLayout(left_panel, stretch=1)

        # ---------------------------------------------------------
        # Right Panel: Status, Progress Bar, and Console / Log Area
        # ---------------------------------------------------------
        right_panel = QVBoxLayout()
        right_panel.setSpacing(8)

        # Status Label and Progress Bar
        status_bar_layout = QHBoxLayout()
        self.status_label = QLabel("Ready")
        status_font = QFont()
        status_font.setBold(True)
        self.status_label.setFont(status_font)
        status_bar_layout.addWidget(self.status_label)
        status_bar_layout.addStretch()

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setFixedHeight(14)
        self.progress_bar.setVisible(False)
        status_bar_layout.addWidget(self.progress_bar)
        right_panel.addLayout(status_bar_layout)

        # Main Log and Output View
        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)
        self.text_edit.setFont(QFont("Consolas", 10))
        right_panel.addWidget(self.text_edit, stretch=1)

        main_layout.addLayout(right_panel, stretch=4)

        # Welcome message in console
        self.append_log("==========================================================================")
        self.append_log(" PHACS-NMR - PHA Polymer Characterization by NMR")
        self.append_log(" Ready. Select an analysis mode or drag and drop a spectrum file.")
        self.append_log("==========================================================================")

    # =========================================================================
    # Logging & Console Management
    # =========================================================================

    def append_log(self, text: str):
        """Append text to the console display safely."""
        self.text_edit.append(text)
        # Auto-scroll to end
        cursor = self.text_edit.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        self.text_edit.setTextCursor(cursor)

    def clear_log(self):
        """Clear the console text area."""
        self.text_edit.clear()

    # =========================================================================
    # Analysis Action Handlers
    # =========================================================================

    def C_NMR_Analysis_click(self):
        """Handle 13C NMR analysis initiation via native file dialog."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select 13C NMR ASCII Spectrum File",
            self.last_directory,
            "ASCII Spectra (*ascii-spec.txt *.txt);;All Files (*)"
        )
        if not file_path:
            self.status_label.setText("13C analysis cancelled.")
            return

        self.last_directory = os.path.dirname(file_path)
        self.prompt_and_run_13c(file_path)

    def prompt_and_run_13c(self, file_path: str):
        """Prompt for noise multiplier and start background 13C analysis."""
        noise_options = [
            "10 (Standard / Noisy signal - Recommended)",
            "5 (Clean signal / High SNR)"
        ]
        choice, ok = QInputDialog.getItem(
            self,
            "13C Noise Level Selection",
            "Select the noise multiplier based on signal quality:",
            noise_options,
            current=0,
            editable=False
        )
        if not ok or not choice:
            self.status_label.setText("13C analysis cancelled.")
            return

        multiple_noise = 10 if "10" in choice else 5
        self.start_worker_task("13C", file_path, {"multiple_noise": multiple_noise})

    def H_NMR_Analysis_click(self):
        """Handle 1H NMR deconvolution analysis initiation via native file dialog."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select 1H NMR Deconvolution Peaks File",
            self.last_directory,
            "Deconvoluted Peak Files (*dconpeaks*.txt *.txt);;All Files (*)"
        )
        if not file_path:
            self.status_label.setText("1H analysis cancelled.")
            return

        self.last_directory = os.path.dirname(file_path)
        self.prompt_and_run_1h(file_path)

    def prompt_and_run_1h(self, file_path: str):
        """Prompt for detection sensitivity and start background 1H analysis."""
        sensitivity, ok = QInputDialog.getDouble(
            self,
            "1H Detection Sensitivity",
            "Monomer detection threshold rate (e.g. 0.01 ignores monomers < 1%):",
            value=0.01,
            min=0.0001,
            max=1.0,
            decimals=4
        )
        if not ok:
            self.status_label.setText("1H analysis cancelled.")
            return

        self.start_worker_task("1H", file_path, {"sensitivity": sensitivity})

    # =========================================================================
    # Drag and Drop Integration
    # =========================================================================

    def is_valid_drag_mime(self, mime_data) -> bool:
        """Check if dropped/dragged mime data contains supported text spectrum files."""
        if mime_data and mime_data.hasUrls():
            urls = mime_data.urls()
            return any(url.toLocalFile().lower().endswith(('.txt', '.csv')) for url in urls)
        return False

    def dragEnterEvent(self, event: QDragEnterEvent):
        """Accept file drop events containing text files."""
        if self.is_valid_drag_mime(event.mimeData()):
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        """Handle dropped spectrum files with smart type detection."""
        mimeData = event.mimeData()
        if not mimeData or not mimeData.hasUrls():
            return

        urls = mimeData.urls()
        dropped_file = urls[0].toLocalFile()
        if not dropped_file or not os.path.isfile(dropped_file):
            return

        self.last_directory = os.path.dirname(dropped_file)
        lower_name = os.path.basename(dropped_file).lower()

        # Keyword auto-detection
        if "dcon" in lower_name or "peak" in lower_name:
            self.append_log(f"\n[Drag & Drop] Detected 1H deconvolution file: {dropped_file}")
            self.prompt_and_run_1h(dropped_file)
        elif "ascii" in lower_name or "spec" in lower_name:
            self.append_log(f"\n[Drag & Drop] Detected 13C ASCII spectrum file: {dropped_file}")
            self.prompt_and_run_13c(dropped_file)
        else:
            # Ambiguous: ask user via clean prompt
            choice, ok = QInputDialog.getItem(
                self,
                "Select Analysis Type",
                f"How would you like to process '{os.path.basename(dropped_file)}'?",
                ["13C NMR Spectrum Analysis", "1H NMR Deconvolution"],
                current=0,
                editable=False
            )
            if not ok or not choice:
                return
            if "13C" in choice:
                self.prompt_and_run_13c(dropped_file)
            else:
                self.prompt_and_run_1h(dropped_file)

    # =========================================================================
    # Worker Thread Lifecycle & Concurrency
    # =========================================================================

    def start_worker_task(self, task_type: str, file_path: str, params: dict):
        """Launch background analysis worker while keeping UI responsive."""
        if self.worker is not None and self.worker.isRunning():
            QMessageBox.warning(
                self,
                "Task Running",
                "An analysis task is already running. Please wait for it to finish."
            )
            return

        # Update UI state for processing
        self._set_ui_busy(True)
        self.btn_open_result_pdf.setVisible(False)
        self.status_label.setText(f"Running {task_type} analysis...")
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate pulsing animation

        # Spawn worker
        self.worker = AnalysisWorker(task_type, file_path, params, parent=self)
        self.worker.progress.connect(self.append_log)
        self.worker.finished.connect(self.on_worker_finished)
        self.worker.failed.connect(self.on_worker_failed)
        self.worker.start()

    def on_worker_finished(self, result: dict):
        """Handle successful analysis completion."""
        self._set_ui_busy(False)
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(100)
        self.status_label.setText("Analysis completed successfully.")

        task_type = result.get("task_type")
        pdf_path = result.get("pdf_path")

        if pdf_path and os.path.isfile(pdf_path):
            self.latest_pdf_path = pdf_path
            self.btn_open_result_pdf.setVisible(True)
            self.append_log(f"\n>>> PDF Report generated: {pdf_path}")
            
            # Offer to open the PDF report
            reply = QMessageBox.information(
                self,
                "Analysis Complete",
                f"{task_type} NMR analysis completed successfully!\n\nReport saved to:\n{os.path.basename(pdf_path)}\n\nWould you like to open the PDF report now?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes
            )
            if reply == QMessageBox.StandardButton.Yes:
                self.open_result_pdf()
        else:
            self.append_log(f"\n>>> {task_type} analysis completed.")

    def on_worker_failed(self, error_msg: str):
        """Handle analysis error."""
        self._set_ui_busy(False)
        self.progress_bar.setVisible(False)
        self.status_label.setText("Analysis failed.")
        QMessageBox.critical(
            self,
            "Analysis Error",
            f"An error occurred during analysis:\n\n{error_msg}\n\nPlease check the spectrum signal quality and file format."
        )

    def _set_ui_busy(self, is_busy: bool):
        """Enable or disable interactive buttons during processing."""
        self.btn_C_NMR_Analysis.setEnabled(not is_busy)
        self.btn_H_NMR_Analysis.setEnabled(not is_busy)
        self.btn_Help.setEnabled(not is_busy)

    # =========================================================================
    # Help, Documentation & PDF Viewing
    # =========================================================================

    def help_clicked(self):
        """Display protocol README text in console and reveal protocol buttons."""
        self.pdf_button.setVisible(True)
        self.button_back.setVisible(True)
        try:
            with open(self.readme_path, "r", encoding="cp1252", errors="replace") as f:
                content = f.read()
            self.text_edit.setPlainText(content)
            self.status_label.setText("Displaying Help & Protocol instructions.")
        except FileNotFoundError:
            QMessageBox.critical(self, "Error", f"Help file not found: {self.readme_path}")

    def back_menu(self):
        """Return from Help view to standard log view without widget layout corruption."""
        self.text_edit.clear()
        self.button_back.setVisible(False)
        self.pdf_button.setVisible(False)
        self.status_label.setText("Ready")
        self.append_log("Returned to main log view.")

    def open_pdf(self):
        """Open informative protocol PDF file in system default viewer."""
        if os.path.isfile(self.pdf_protocol_path):
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.pdf_protocol_path))
        else:
            QMessageBox.warning(self, "Warning", f"Protocol PDF not found at {self.pdf_protocol_path}")

    def open_result_pdf(self):
        """Open the latest generated analysis PDF in system default viewer."""
        if self.latest_pdf_path and os.path.isfile(self.latest_pdf_path):
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.latest_pdf_path))
        else:
            QMessageBox.warning(self, "Warning", "No PDF report is currently available to open.")

    # =========================================================================
    # Window Close Event & Cleanup
    # =========================================================================

    def closeEvent(self, event):
        """Confirm exit if worker is busy and restore sys.stdout."""
        if self.worker is not None and self.worker.isRunning():
            reply = QMessageBox.question(
                self,
                "Analysis In Progress",
                "An analysis is currently running in the background. Are you sure you want to abort and exit?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            if reply != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.worker.terminate()
            self.worker.wait(1000)

        # Restore original stdout and stderr
        if sys.stdout is self.stdout_redirector:
            sys.stdout = self.original_stdout
        if sys.stderr is self.stderr_redirector:
            sys.stderr = self.original_stderr

        event.accept()


# =============================================================================
# Application Entry Point
# =============================================================================

def main():
    """Application entry point for PHACS-NMR GUI."""
    # Ensure Windows taskbar groups and displays the application icon properly
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("bioeconomyscienceinstitute.phacsnmr.app.2.0")
        except Exception:
            pass

    app = cast(QApplication, QApplication.instance())
    if app is None:
        app = QApplication(sys.argv)
    app.setApplicationName("PHACS-NMR")

    icon_path = _find_resource("phacs_nmr.ico")
    if os.path.isfile(icon_path):
        app.setWindowIcon(QIcon(icon_path))

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
