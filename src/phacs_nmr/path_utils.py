"""
Path utilities for PHACS-NMR.
Provides robust resolution of sample names and output filenames across
Bruker TopSpin folder hierarchies, nested test directories, and flat files.
"""
from pathlib import Path


def is_drive_or_root(part: str) -> bool:
    """Returns True if a path part represents a root or drive specification."""
    if not part:
        return True
    clean = part.rstrip('\\/ ')
    return clean.endswith(':') or clean in ('', '/', '\\')


DEFAULT_TOPSPIN_STEMS = {"ascii-spec", "dconpeaks"}


def resolve_sample_name(filepath: str | Path) -> str:
    """
    Extracts the sample name from a file path:
    - If 'pdata' is present in the path parts and len(parts) >= 5:
      returns parts[-5] (the Bruker TopSpin sample directory).
    - Otherwise (if 'pdata' is not present or path has < 5 parts / drive root):
      returns the file stem.
    Guarantees that drive letters ('D:', 'D:\\', 'C:\\') are never returned as sample names.
    """
    p = Path(filepath).resolve()
    parts = p.parts
    lower_parts = [part.lower() for part in parts]

    # 1. TopSpin pdata hierarchy: <sample>/<expno>/pdata/<procno>/<file>
    if 'pdata' in lower_parts:
        if len(parts) >= 5 and not is_drive_or_root(parts[-5]):
            return parts[-5]

    # 2. No 'pdata' or root fallback: use file stem
    return p.stem


def resolve_output_prefix(filepath: str | Path) -> str:
    """
    Resolves the filename prefix for generated results and reports:
    - If 'pdata' in path:
        - sample_name if stem in DEFAULT_TOPSPIN_STEMS ('ascii-spec', 'dconpeaks')
        - f"{sample_name}_{p.stem}" if custom non-default stem
    - If 'pdata' not in path:
        - file stem
    """
    p = Path(filepath).resolve()
    sample_name = resolve_sample_name(p)
    lower_parts = [part.lower() for part in p.parts]

    if 'pdata' in lower_parts:
        if sample_name == p.stem or p.stem.lower() in DEFAULT_TOPSPIN_STEMS:
            return sample_name
        return f"{sample_name}_{p.stem}"
    return p.stem
