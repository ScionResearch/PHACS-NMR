from pathlib import Path
import numpy as np


def generate_synthetic_phbhhxho_mcl_ascii(target_file: Path, n_points: int = 65536):
    """
    Generate a high-fidelity synthetic 13C ASCII spectrum for PHBHHxHO_mcl
    containing the characteristic peaks and region resonances needed by table1.
    """
    ppm = np.linspace(220.0, -20.0, n_points)
    # Baseline noise ~ 1e5
    np.random.seed(42)
    intensity = np.random.normal(0, 1e5, n_points)

    def add_peak(center_ppm, height=1e8, fwhm_ppm=0.02):
        gamma = fwhm_ppm / 2.0
        return height * (gamma**2 / ((ppm - center_ppm)**2 + gamma**2))

    # 1. Characteristic peaks for find_polymer_ID: [0, 1, 0, 1, 1, 1]
    intensity += add_peak(19.75, height=2e8)   # B4 (19.5-20.0) -> B=1
    intensity += add_peak(13.80, height=2e8)   # X6 (13.7-13.9) -> X=1
    intensity += add_peak(22.55, height=2e8)   # O7 (22.4-22.7) -> O=1
    intensity += add_peak(25.10, height=1.5e8) # mcl (25.0-25.2) -> mcl=1

    # 2. Resonances centered within each table1 region:
    intensity += add_peak(169.30, height=3e8)  # B1X1O1 (169.0-169.6)
    intensity += add_peak(70.85, height=1.5e8) # O3 (70.7-71.0)
    intensity += add_peak(70.60, height=1.5e8) # X3 (70.5-70.7)
    intensity += add_peak(67.60, height=2e8)   # B3 (67.4-67.8)
    intensity += add_peak(40.85, height=2e8)   # B2 (40.0-42.0)
    intensity += add_peak(39.10, height=1.5e8) # O2X2 (38.5-39.5)
    intensity += add_peak(36.00, height=1.5e8) # X4 (35.5-36.5)
    intensity += add_peak(33.75, height=1.5e8) # O4 (33.65-33.85)
    intensity += add_peak(31.50, height=1.5e8) # O5 (30.6-31.8, kept away from F2 30.5-30.9)
    intensity += add_peak(24.75, height=1.5e8) # O6 (24.6-24.9)
    intensity += add_peak(18.40, height=1.5e8) # X5 (18.2-18.6)
    intensity += add_peak(14.12, height=1.2e8) # mcl (14.093-14.16)
    intensity += add_peak(14.00, height=1.5e8) # O8 (13.9-14.1)

    # Format as standard 4-column TopSpin convbin2asc export
    lines = []
    for idx in range(n_points):
        # 1, intensity, 0.0, ppm
        lines.append(f"{idx + 1}, {intensity[idx]:.1f}, 0.0, {ppm[idx]:.6f}\n")

    target_file.parent.mkdir(parents=True, exist_ok=True)
    target_file.write_text("".join(lines), encoding="cp1252")
    return target_file
