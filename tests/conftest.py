"""
Pytest configuration and session-wide fixtures for PHACS-NMR.

Provides:
- .env loading for external historical benchmark paths
- sys.path configuration ensuring project root is prioritized
- benchmark_data_dir: skips tests if local external benchmark disk is absent
- global_software_dir: skips tests if external 1H deconvolution benchmark disk is absent
- synthetic_spectrum: in-memory synthetic 13C spectrum array with known Lorentzian peaks
- cleanup_matplotlib_figures: autouse fixture closing all matplotlib figures after each test
- qapp: session-scoped offscreen QApplication instance for headless PyQt6 GUI testing
"""
import os
import sys
from pathlib import Path

# Load local .env file if present (ignoring comments and empty lines)
_env_path = Path(__file__).resolve().parent.parent / ".env"
if _env_path.is_file():
    for _line in _env_path.read_text(encoding="utf-8").splitlines():
        _line = _line.strip()
        if _line and not _line.startswith("#") and "=" in _line:
            _k, _v = _line.split("=", 1)
            os.environ.setdefault(_k.strip(), _v.strip().strip('"').strip("'"))

# Ensure src directory is on sys.path before third-party imports
SRC_DIR = Path(__file__).resolve().parent.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pytest
import numpy as np

# Resolve external historical benchmark directories from environment variables
_bench_env = os.environ.get("PHACS_BENCHMARK_DIR")
BENCHMARK_DATA_DIR = Path(_bench_env).resolve() if _bench_env else None

_global_env = os.environ.get("PHACS_GLOBAL_SOFTWARE_DIR")
GLOBAL_SOFTWARE_DIR = Path(_global_env).resolve() if _global_env else None


@pytest.fixture
def benchmark_data_dir():
    """Return Path to benchmark directory if it exists, else skip."""
    if BENCHMARK_DATA_DIR is None or not BENCHMARK_DATA_DIR.exists():
        pytest.skip(f"Benchmark directory not configured or found: {BENCHMARK_DATA_DIR}")
    return BENCHMARK_DATA_DIR


@pytest.fixture
def global_software_dir():
    """Return Path to global software additional info directory if it exists, else skip."""
    if GLOBAL_SOFTWARE_DIR is None or not GLOBAL_SOFTWARE_DIR.exists():
        pytest.skip(f"Global software directory not configured or found: {GLOBAL_SOFTWARE_DIR}")
    return GLOBAL_SOFTWARE_DIR


@pytest.fixture
def synthetic_spectrum():
    """Generate a simple synthetic 13C spectrum array with baseline and known peaks."""
    # PPM scale from 220 down to -20 with 65536 points
    n_points = 65536
    ppm_scale = np.linspace(220.0, -20.0, n_points)
    intensities = np.random.normal(0, 1e6, n_points)

    def add_lorentzian(center_ppm, height, fwhm_ppm):
        gamma = fwhm_ppm / 2.0
        return height * (gamma**2 / ((ppm_scale - center_ppm)**2 + gamma**2))

    # Add characteristic peaks for PHBV:
    # B4: 19.74 ppm
    intensities += add_lorentzian(19.74, 5e7, 0.05)
    # V5: 9.32 ppm
    intensities += add_lorentzian(9.32, 3e7, 0.05)
    # B1: 169.14 ppm
    intensities += add_lorentzian(169.14, 4e7, 0.05)
    # V1: 169.51 ppm
    intensities += add_lorentzian(169.51, 2e7, 0.05)

    return ppm_scale, intensities


@pytest.fixture(autouse=True)
def cleanup_matplotlib_figures():
    """Ensure all matplotlib figures are closed after each test to prevent memory accumulation."""
    import matplotlib.pyplot as plt
    yield
    plt.close('all')


@pytest.fixture(scope="session")
def qapp():
    """Create or retrieve a headless offscreen QApplication instance for tests."""
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        app = QApplication(["--platform", "offscreen"])
    yield app


