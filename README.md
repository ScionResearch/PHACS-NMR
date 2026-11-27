# PHACS-NMR

**PHACS-NMR - PHA Polymer Characterization by NMR**

`PHACS-NMR` is a scientific software package developed at the Bioeconomy Science Institute (New Zealand Institute for Bioeconomy Science Limited) for the automated compositional and microstructural characterization of bacterially synthesized polyhydroxyalkanoate (PHA) biopolymers from Bruker TopSpin NMR data.

> [!NOTE]
> The methodology, peak deconvolution models, and validation are described in: 
> Fornes, C., Pion, F., Thomas, H., Schou, W., et al. Automated Proton and Carbon NMR Analysis of Polyhydroxyalkanoate Composition and Sequence Distribution. *Submitted to Analytical Chemistry*.

---

## Key Capabilities

1. **$^1\text{H}$ NMR Fast Monomer Composition Screening (`dconpeaks.txt`)**:
   - Quantifies monomer mole fractions in $< 10\text{ minutes}$ (including acquisition) from Bruker TopSpin deconvolution peak tables (`dcon`).
   - Supports 9 monomer classes: **4HB** (`F`), **3HB** (`B`), **3HV** (`V`), **3HHx** (`X`), **3HHp** (`Hp`), **3HO** (`O`), **3HN** (`N`), and long-chain **MCL** ($>\text{C}_9$) units.
   - Applies empirically calibrated multiplet and stoichiometric correction factors and exports a timestamped CSV report.

2. **$^{13}\text{C}$ NMR Quantitative Copolymer Microstructure Analysis (`ascii-spec.txt`)**:
   - Automatically detects the PHA copolymer family from diagnostic $^{13}\text{C}$ chemical shifts (`PHBV`, `P34HB`, `PHBHHx`, `PHBHHxHO_MCL`, or `PHB` homopolymer).
   - Performs bounded Levenberg-Marquardt Lorentzian lineshape fitting (`nmrglue`) across all carbonyl, methine, methylene, and methyl multiplets.
   - Computes monad ($F_X$), diad ($F_{XY}$), and triad ($F_{XYZ}$) sequence distributions across homologous carbon pairs ($C_1\text{--}C_1$, $C_2\text{--}C_2$, $C_3\text{--}C_3$, $C_4\text{--}C_4$).
   - Evaluates copolymer randomness and blockiness metrics:
     - **Kamiya $D$ Parameter** ($D = 1$ random, $D > 1$ blocky/blend, $D < 1$ alternating)
     - **Stanford ChemE 160 $\chi$ Parameter** ($\chi = 1$ random, $\chi < 1$ blocky, $\chi > 1$ alternating)
     - **Number-Average Block Lengths** ($\bar{l}_A, \bar{l}_B$)
   - Fits three statistical copolymerization models: **Model 1** (Bernoullian random), **Model 2** (1st-order Markovian), and **Model 3** (Two-component Bernoullian blend via `scipy.optimize.fsolve`).
   - Generates a 5-figure multi-page **PDF graphical report** and a comprehensive **CSV summary table**.

---

## Documentation

- 📘 **[Copolymer Microstructure & Statistical Models](docs/copolymer_models.md)**: Mathematical derivations, KaTeX equations, physical interpretations of $D$, $\chi$, $\bar{l}_A, \bar{l}_B$, and Models 1–3 grounded in peer-reviewed literature.
- 📗 **[Chemical Shift Assignments, Calibration & Spectral Processing](docs/calibration_and_assignments.md)**: Complete $^1\text{H}$ and $^{13}\text{C}$ chemical shift tables, empirical deconvolution calibration constants ($0.5437$, $0.8401$, etc.), and TopSpin acquisition/preprocessing protocols.
- 📄 **[Laboratory Software Protocol v2 (PDF)](Software%20Protocol%20v2.pdf)**: Step-by-step wet-lab sample preparation, Bruker Avance NEO 600 MHz acquisition parameters (`zg30` / `zgig30`), and TopSpin 4.5.0 processing screenshots.

---

## Installation & Setup

### Prerequisites
- **Python**: `>=3.10, <3.12` (tested on Python 3.11)
- **Operating System**: Windows, macOS, or Linux

### Option A: Using `uv` (Recommended)
```pwsh
# From the repository root
uv venv .venv --python 3.11
.\.venv\Scripts\activate
uv pip install -e ".[dev]"
```

### Option B: Using Standard `pip`
```pwsh
# From the repository root
python -m venv .venv
.\.venv\Scripts\activate
pip install -e ".[dev]"
```

---

## Building & Packaging

To create redistributable Python wheel (`.whl`) and source distribution (`.tar.gz`) packages:

### Option A: Using `uv` (Recommended)
```pwsh
uv build
```

### Option B: Using Standard `build`
```pwsh
pip install build
python -m build
```

Both commands package the project defined in `pyproject.toml` and output distribution files into the `dist/` folder:
- `dist/phacs_nmr-2.0.0-py3-none-any.whl` (Binary wheel)
- `dist/phacs_nmr-2.0.0.tar.gz` (Source distribution)

### Installing the Built Wheel
The resulting wheel can be installed into any compatible Python environment without needing the source repository:
```pwsh
pip install dist/phacs_nmr-2.0.0-py3-none-any.whl
```

---

## Usage Guide

### 1. Graphical User Interface (GUI)

Launch the desktop application using any of the following commands:
```pwsh
# Installed console script (recommended):
phacs-nmr-gui

# Or via Python module execution:
python -m phacs_nmr.gui
# (or simply: python -m phacs_nmr)
```

#### GUI Workflows:
- **Drag & Drop**: Drag any `.txt` file (`ascii-spec.txt` or `dconpeaks.txt`) directly onto the window. The application automatically detects whether the file is a $^{13}\text{C}$ ASCII spectrum or a $^1\text{H}$ deconvolution peak list based on filename keywords and prompts only if ambiguous.
- **C NMR Analysis (`ascii-spec.txt`)**:
  1. Click **C NMR Analysis** and select your TopSpin `ascii-spec.txt` file via the native file picker.
  2. Select the noise multiplier:
     - **`10 - Noisy signal (standard)`**: Recommended for routine $2048$-scan spectra or samples with minor baseline impurities.
     - **`5 - Clean signal (high S/N)`**: Recommended for high-scan spectra ($4096\text{--}16384$ scans) or when quantifying minor monomer fractions ($< 3\%$).
  3. The analysis runs in a non-blocking background thread (`QThread`) with live progress streamed to the console pane.
  4. Upon completion, click **Yes** (or click the **Open Latest PDF Report** button at any time) to view the generated 5-figure PDF report (`{sample}_13C_n{noise}_Software_results_{timestamp}.pdf`).
- **H NMR Analysis (`dconpeaks.txt`)**:
  1. Click **H NMR Analysis** and select your TopSpin `dconpeaks.txt` file.
  2. Confirm or adjust the detection sensitivity threshold (default `0.01`, i.e., $1\%$).
  3. View the formatted monomer composition table directly in the log window; the CSV report (`{sample}_H_s{sensitivity}_Software_results_{timestamp}.csv`) is saved automatically next to the input file.

---

### 2. Command-Line Processing (CLI)

`phacs-nmr` operates as a unified command-line tool that automatically handles single spectrum files or batch processing across entire directories.

#### A. Single Spectrum Analysis
To run polymer auto-detection, deconvolution, and reporting on a single $^{13}\text{C}$ ASCII spectrum file:

```pwsh
phacs-nmr <path_to_txt_file> [noise_level]
# Or using the optional flag:
phacs-nmr <path_to_txt_file> -n <5|10>
# Or via python -m:
python -m phacs_nmr <path_to_txt_file> [noise_level]
```

**Arguments**:
- `<path_to_txt_file>`: Path to the input NMR ASCII text file (e.g. `ascii-spec.txt`).
- `[noise_level]` or `-n/--noise`: `5` (clean) or `10` (noisy, default).

#### B. Batch Directory Processing
To process an entire directory of $^{13}\text{C}$ `*.txt` spectrum files headlessly without launching the GUI:

```pwsh
phacs-nmr <path_to_folder> [noise_level]
```

**Arguments**:
- `<path_to_folder>`: Directory containing `.txt` files exported via TopSpin `convbin2asc`.
- `[noise_level]`: `5` (clean spectra) or `10` (noisy spectra, default).

**Example**:
```pwsh
phacs-nmr "D:\NMR_Data\Batch_2026_10" 5
```
Each spectrum is processed independently with automatic error isolation; generated `.pdf` and `.csv` reports are written alongside each input `.txt` file.

---

### 3. Programmatic Python API

You can also invoke the analysis engines directly from Python scripts or Jupyter notebooks:

```python
# 13C NMR Copolymer Microstructure Analysis
import phacs_nmr

figures = phacs_nmr.start_analysis(
    input_file="path/to/sample/11/pdata/1/ascii-spec.txt",
    multiple_noise=5,
)

# 1H NMR Deconvolution Composition Analysis
from phacs_nmr.core.io import results_dconv

headers, p_id, polymer_type, areas, intensities = results_dconv(
    filename="path/to/sample/10/pdata/1/dconpeaks.txt",
    sensitivity=0.01,
)
```

---

## Project Architecture

```text
phacs-nmr/
├── src/
│   └── phacs_nmr/                      # Primary application package
│       ├── __init__.py                 # Public package interface & version
│       ├── __main__.py                 # `python -m phacs_nmr` entrypoint
│       ├── cli.py                      # CLI runner & headless batch processor
│       ├── gui.py                      # Modernized PyQt6 GUI with QThread worker
│       ├── path_utils.py               # TopSpin path & sample name resolution
│       ├── assets/                     # Application icons (.ico, .png, .svg)
│       └── core/                       # Scientific processing engine
│           ├── __init__.py
│           ├── config.py               # Declarative PolymerProfile definitions & ppm limits
│           ├── engine.py               # End-to-end 13C analysis, plotting, PDF & CSV export
│           ├── io.py                   # TopSpin ascii-spec.txt & dconpeaks.txt parsers & 1H calibration
│           ├── models.py               # Kamiya Models 1, 2, 3, D parameter, Stanford chi & block lengths
│           ├── noise.py                # RMS baseline noise estimation & ppm index scaling
│           └── peak_fitting.py         # Spectral region slicing & nmrglue Lorentzian fitting
├── docs/                               # Scientific & technical documentation
│   ├── copolymer_models.md             # Mathematical derivations & literature mapping
│   └── calibration_and_assignments.md     # Chemical shifts, calibration factors & TopSpin protocol
├── tests/                              # Automated pytest suite (79 unit & integration tests)
│   ├── conftest.py                     # Shared fixtures & headless PyQt6 QApplication
│   ├── fixtures/                       # Synthetic spectrum generators
│   ├── integration/                    # Golden master parity, batch runner & pipeline tests
│   └── unit/                           # Unit tests for models, I/O, noise, paths, GUI & QThread worker
├── Software Protocol v2.pdf            # Laboratory NMR acquisition & TopSpin user protocol
├── pyproject.toml                      # Project dependencies, build config & pytest settings
└── README.md
```

---

## Testing & Verification

The repository includes a comprehensive `pytest` suite verifying unit logic, GUI concurrency, and bit-for-bit numerical parity against golden reference outputs:

```pwsh
# Run the full test suite
.\.venv\Scripts\pytest.exe -v

# Run with statement coverage report
.\.venv\Scripts\pytest.exe --cov=core --cov=gui --cov-report=term-missing
```

---

## Credits & References

### Institution & Contributors
- **Institution**: Bioeconomy Science Institute (New Zealand Institute for Bioeconomy Science Limited).
- **Authors & Contributors**: Clément Fornes, Samuel Damesin, Ralf Gommers, Wayne Schou.
- **Scientific Advisors**: Mathias Sorieul, Christophe Collet, Andrew Vogt, Marion Sanglard, Marc Gaugler, Stefan Hill.

### Associated Publication
- Fornes, C., Pion, F., Thomas, H., Schou, W., et al. Automated Proton and Carbon NMR Analysis of Polyhydroxyalkanoate Composition and Sequence Distribution. *Submitted to Analytical Chemistry*.

### Key Literature
1. Kamiya, N.; Yamamoto, Y.; Inoue, Y.; Chûjô, R.; Doi, Y. Microstructure of Bacterially Synthesized Poly(3-hydroxybutyrate-co-3-hydroxyvalerate). *Macromolecules* **1989**, *22* (4), 1676–1682. DOI: [10.1021/ma00194a030](https://doi.org/10.1021/ma00194a030).
2. Ivanova, G.; Serafim, L. S.; Lemos, P. C.; Ramos, A. M.; Reis, M. A. M.; Cabrita, E. J. Influence of feeding strategies of mixed microbial cultures on the chemical composition and microstructure of copolyesters P(3HB-co-3HV) analyzed by NMR and statistical analysis. *Magn. Reson. Chem.* **2009**, *47* (6), 497–504. DOI: [10.1002/mrc.2423](https://doi.org/10.1002/mrc.2423).
3. Frank, C. W. *Lecture 7: Statistics of Chain Copolymerization*; Chemical Engineering 160/260 (Polymer Science and Engineering), Stanford University, **2001**.
4. Shapiro, Yu. E. Analysis of Chain Microstructure by $^1\text{H}$ and $^{13}\text{C}$ NMR Spectroscopy. *Bull. Magn. Reson.* **1985**, *7* (1), 27–58.
5. Suganuma, K.; Asakura, T.; Oshimura, M.; Hirano, T.; Ute, K.; Cheng, H. N. NMR Analysis of Poly(Lactic Acid) via Statistical Models. *Polymers* **2019**, *11* (4), 725. DOI: [10.3390/polym11040725](https://doi.org/10.3390/polym11040725).
6. Fornes, C. *Protocol for using PHA software: $^1\text{H}$ NMR characterization of monomer composition and $^{13}\text{C}$ NMR characterization of copolymer microstructure*; Technical Protocol v2, December **2025**.

---

## Third-Party Package Acknowledgments & Citations

`PHACS-NMR` is built on top of the scientific Python ecosystem. In accordance with the citation requests of their respective authors and development teams, the third-party packages utilized by `PHACS-NMR` are acknowledged and cited below:

- **`nmrglue`** (BSD 3-Clause): J. J. Helmus, C. P. Jaroniec, *nmrglue: An open source Python package for the analysis of multidimensional NMR data*, *Journal of Biomolecular NMR* **2013**, *55* (4), 355–367. DOI: [10.1007/s10858-013-9718-x](https://doi.org/10.1007/s10858-013-9718-x).
- **`SciPy`** (BSD 3-Clause): P. Virtanen, R. Gommers, T. E. Oliphant, M. Haberland, T. Reddy, D. Cournapeau, E. Burovski, P. Peterson, W. Weckesser, J. Bright, S. J. van der Walt, M. Brett, J. Wilson, K. J. Millman, N. Mayorov, A. R. J. Nelson, E. Jones, R. Kern, E. Larson, C. J. Carey, İ. Polat, Y. Feng, E. W. Moore, J. VanderPlas, D. Laxalde, J. Perktold, R. Cimrman, I. Henriksen, E. A. Quintero, C. R. Harris, A. M. Archibald, A. H. Ribeiro, F. Pedregosa, P. van Mulbregt, and SciPy 1.0 Contributors. *SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python*, *Nature Methods* **2020**, *17* (3), 261–272. DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).
- **`NumPy`** (BSD 3-Clause): C. R. Harris, K. J. Millman, S. J. van der Walt, R. Gommers, P. Virtanen, D. Cournapeau, E. Wieser, J. Taylor, S. Berg, N. J. Smith, R. Kern, M. Picus, S. Hoyer, M. H. van Kerkwijk, M. Brett, A. Haldane, J. Fernández del Río, M. Wiebe, P. Peterson, P. Gérard-Marchant, K. Sheppard, T. Reddy, W. Weckesser, H. Abbasi, C. Gohlke, and T. E. Oliphant. *Array programming with NumPy*, *Nature* **2020**, *585*, 357–362. DOI: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2).
- **`Matplotlib`** (PSF / Matplotlib): J. D. Hunter. *Matplotlib: A 2D Graphics Environment*, *Computing in Science & Engineering* **2007**, *9* (3), 90–95. DOI: [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55).
- **`pandas`** (BSD 3-Clause): The pandas development team. *pandas-dev/pandas: Pandas*, Zenodo, **2020**. DOI: [10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134); and W. McKinney, *Data Structures for Statistical Computing in Python*, *Proceedings of the 9th Python in Science Conference* **2010**, 51–56. DOI: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a).
- **`seaborn`** (BSD 3-Clause): M. L. Waskom. *seaborn: statistical data visualization*, *Journal of Open Source Software* **2021**, *6* (60), 3021. DOI: [10.21105/joss.03021](https://doi.org/10.21105/joss.03021).
- **`ReportLab`** (BSD License): ReportLab Inc. *ReportLab PDF Library*, Wimbledon, London, UK. URL: [https://www.reportlab.com/](https://www.reportlab.com/).
- **`PyQt6`** (GNU GPL v3 / Commercial): Riverbank Computing Limited. *PyQt6: Python Bindings for Qt 6*. URL: [https://www.riverbankcomputing.com/software/pyqt/](https://www.riverbankcomputing.com/software/pyqt/).
- **`tabulate`** (MIT License): S. Astanin. *tabulate: Pretty-print tabular data in Python*. URL: [https://github.com/astanin/python-tabulate](https://github.com/astanin/python-tabulate).
- **`pytest` & `pytest-cov`** (MIT / BSD): H. Krekel et al. *pytest: Simple powerful testing with Python*. URL: [https://github.com/pytest-dev/pytest](https://github.com/pytest-dev/pytest).

---

## License and Disclaimer  
PHACS-NMR is provided as research software. Users are responsible for assessing its suitability for their particular purposes and for independently validating results where appropriate. The software is provided subject to the warranty and liability provisions of the GNU GPLv3.

See the full license text in the [`LICENSE`](LICENSE) file.