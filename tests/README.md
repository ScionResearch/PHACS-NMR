# PHACS-NMR Test Suite Documentation

This directory contains the automated test suite for **PHACS-NMR**, an automated NMR analysis and copolymer microstructure quantification platform for polyhydroxyalkanoates (PHAs).

---

## 1. Directory Structure

```
tests/
├── README.md                      # This documentation
├── conftest.py                    # Global fixtures (QApplication, Matplotlib cleanup, synthetic spectra)
├── fixtures/                      # Reusable synthetic spectrum and peak generators
│   └── synthetic_mcl_spectra.py   # High-fidelity synthetic 13C ASCII spectrum for MCL terpolymers
├── integration/                   # End-to-end and regression tests (historical golden parity & pipelines)
│   ├── test_batch_runner.py       # Batch processing execution and error isolation
│   ├── test_dconv_golden.py       # 1H TopSpin deconvolution parity against 2025 CSV baselines
│   ├── test_golden_parity.py      # 13C lineshape & Kamiya parity against 2026 CSV baselines
│   └── test_mcl_pipeline.py       # Full synthetic MCL terpolymer analysis pipeline
└── unit/                          # Fast, self-contained unit and invariant tests (< 6s)
    ├── test_config_profiles.py    # Declarative PolymerProfiles and chemical shift limit builders
    ├── test_copolymer_models.py   # Kamiya Models 1-3, D, chi, block lengths, and physical invariants
    ├── test_core_engine_ratios.py # Peak integration, area normalization, and CSV export
    ├── test_core_io_robustness.py # ASCII/dcon file parsers, corrupt line guards, sensitivity filtering
    ├── test_dconv_all_monomers.py # 1H deconvolution monomer classification for terpolymers & MCL
    ├── test_dconv_correction.py   # Monomer limits intervals and PHBV peak detection
    ├── test_gui_window.py         # PyQt6 MainWindow layout, drag-and-drop, and log consoles
    ├── test_gui_worker.py         # PyQt6 AnalysisWorker background threading and signal lifecycle
    ├── test_noise_and_peaks.py    # Baseline noise estimation (-5 to -15 ppm) and ppm scale indexing
    ├── test_path_utils.py         # Bruker TopSpin directory hierarchy and sample name resolution
    └── test_synthetic_deconvolution.py # Lineshape fitting on synthetic Lorentzian peaks & noise rejection
```

---

## 2. Test Execution Quick Start

All tests are managed via `uv` or standard `pytest`.

### Running Fast Unit Tests (< 6 seconds)
These tests do **not** require any local benchmark files or external data drives:
```bash
uv run pytest tests/unit/ -v
```

### Running with Code Coverage
To generate a terminal coverage report for the `core/` analytical engine:
```bash
uv run coverage run -m pytest tests/unit/
uv run coverage report --include="src/phacs_nmr/core/*"
```

### Running Full Test Suite (Unit + Integration, ~60 seconds)
```bash
uv run pytest
```
> [!NOTE]
> Tests in `tests/integration/test_golden_parity.py` and `tests/integration/test_dconv_golden.py` compare newly generated analysis CSVs against historical benchmark CSVs configured via environment variables (`PHACS_BENCHMARK_DIR` and `PHACS_GLOBAL_SOFTWARE_DIR` in `.env`). If these environment variables or directories are not present, those specific golden tests are skipped automatically via `@pytest.mark.skipif`. Copy `.env.example` to `.env` to configure paths on your local machine.

---

## 3. Testing Methodologies & Categories

The suite employs two complementary testing paradigms:

### A. Mathematical Invariant & Property Tests (`tests/unit/`)
Rather than relying on legacy software output, invariant tests verify fundamental physical and statistical laws of polymer chemistry:
1. **Monomer Conservation:** Monomer fractions must sum to unity ($\sum f_i = 1.0$).
2. **Diad Normalization:** $F_{AA} + F_{AB} + F_{BA} + F_{BB} = 1.0$.
3. **Triad Marginals:** Triad probabilities must sum to their parent diads ($F_{AAA} + F_{AAB} = F_{AA}$).
4. **Bernoullian Randomness Invariants:** For any ideal random distribution, Kamiya $D \equiv 1.0$ and Stanford $\chi \equiv 1.0$ within $10^{-10}$ floating-point precision.
5. **Model 3 Blend Parameter Recovery:** Synthetic triad mixtures with known true copolymer compositions ($A, B$) and blend fraction ($X$) must be recovered accurately by the `scipy.optimize.fsolve` solver.
6. **Numerical Guards & Sentinels:** Homopolymers with zero cross-diads cleanly return $D = 100.0$ and $\chi = \text{NaN}$ without raising `ZeroDivisionError`.

### B. Synthetic Deconvolution & Ground-Truth Recovery (`tests/unit/test_synthetic_deconvolution.py`)
* Synthesizes in-memory NMR spectra with known Lorentzian/Voigt peaks (analytical centers, widths, heights, and areas) combined with Gaussian noise.
* Uses the `preloaded_xy` interface of `split_graphic_regions` to bypass disk I/O and verify peak picking and curve fitting recovery against analytical ground truth.
* Tests baseline noise rejection ($10\sigma$) and absent monomer handling ($r == 0$).

### C. Golden Parity Regression Tests (`tests/integration/`)
* Runs full analysis on raw TopSpin spectrometer files and asserts parity against historical 2025 and 2026 CSV outputs.
* Ensures strict backwards compatibility during refactoring.

---

## 4. Scientific Literature Citations

Test calculations, parameter definitions, and expected values are derived from:

1. **Kamiya, N., et al. (1989):**
   * *Microstructure of bacterially synthesized poly(3-hydroxybutyrate-co-3-hydroxyvalerate).* Macromolecules, 22(4), 1676-1682.
   * Model 1 (Bernoullian random copolymer, Eqs 1 & 3).
   * Model 2 (First-order Markovian copolymer, Eq 4).
   * Model 3 (Copolymer blend, Eq 6).
   * Kamiya Parameter: $D = \frac{F_{AA} F_{BB}}{F_{AB} F_{BA}}$ (Eq 2).
2. **Stanford ChemE 160 Lecture Notes:**
   * Sequence distribution and block length statistics:
     $$\chi = \frac{P_2\{BA\}}{P_1\{A\} P_1\{B\}}, \quad \bar{L}_A = \frac{A^*}{A^*B}, \quad \bar{L}_B = \frac{B^*}{B^*A}$$
3. **Bruker TopSpin convbin2asc & Deconvolution Protocol:**
   * TopSpin 4-column ASCII export (`convbin2asc`): `[index, intensity, 0.0, ppm]`.
   * TopSpin `dconpeaks.txt` peak report: 8 header lines followed by `[ppm, 0, 0, 0, intensity, area, fwhm]`.
4. **Core Scientific Packages:**
   * **nmrglue**: Helmus & Jaroniec, *J. Biomol. NMR* 2013, 55, 355–367. DOI: 10.1007/s10858-013-9718-x.
   * **SciPy**: Virtanen et al., *Nature Methods* 2020, 17, 261–272. DOI: 10.1038/s41592-019-0686-2.
   * **NumPy**: Harris et al., *Nature* 2020, 585, 357–362. DOI: 10.1038/s41586-020-2649-2.
   * **Matplotlib**: Hunter, *Comput. Sci. Eng.* 2007, 9, 90–95. DOI: 10.1109/MCSE.2007.55.
   * **pandas**: McKinney, *Proc. 9th Python Sci. Conf.* 2010, 51–56; pandas dev team, Zenodo 2020.
   * **seaborn**: Waskom, *J. Open Source Softw.* 2021, 6(60), 3021. DOI: 10.21105/joss.03021.

---

## 5. Guidelines for Writing New Tests

When adding new tests to this repository, please adhere to the following standards:

1. **Prefer Synthetic In-Memory Testing:**
   * Do not rely on external disk files for unit tests. Use programmatic arrays, `preloaded_xy`, or `tmp_path`.
2. **Follow the Given / When / Then Pattern:**
   * **Given:** Document the input state (e.g. chemical shifts, noise sigma, monomer fractions).
   * **When:** Call the specific function under test.
   * **Then:** Assert expected results with clear domain failure messages.
3. **Provide Informative Assertion Messages:**
   * Always include context in assertion messages explaining *why* the assertion matters:
     ```python
     assert np.isclose(d_val, 1.0, atol=1e-6), (
         f"Kamiya D must equal 1.0 for ideal random copolymer, got {d_val}"
     )
     ```
4. **Isolate GUI and Matplotlib State:**
   * Use the `qapp` session fixture for any PyQt6 widget tests.
   * Matplotlib figures are automatically closed after each test via the `cleanup_matplotlib_figures` fixture in `conftest.py`.
