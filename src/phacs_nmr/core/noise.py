"""
Noise estimation and ppm scale indexing utilities for PHACS-NMR.

References:
    Harris, C.R. et al. (2020). Array programming with NumPy. Nature, 585, 357-362.
"""
import numpy as np


def inv_scale(array: np.ndarray | list, value: float) -> int:
    """Find closest index in array to the given value."""
    arr = np.asarray(array)
    return int(np.abs(np.subtract.outer(arr, value)).argmin(0))


def calc_noise_level(
    ppm_scale: np.ndarray | list,
    data: np.ndarray | list,
    noise_range: tuple[float, float] = (-5.0, -15.0),
    verbose: bool = True,
) -> float:
    """
    Calculate signal noise using the average absolute value of the region
    (e.g., -5 to -15 ppm) where no NMR signal should ever appear.
    """
    ppm = np.asarray(ppm_scale)
    y = np.asarray(data)

    idx_start = inv_scale(ppm, noise_range[0])
    idx_end = inv_scale(ppm, noise_range[1])
    slice_start = min(idx_start, idx_end)
    slice_end = max(idx_start, idx_end)

    noise_slice = y[slice_start:slice_end]
    noise = float(np.mean(np.abs(noise_slice)))
    noise_rms = float(np.sqrt(np.mean(np.square(noise_slice))))

    if verbose:
        mean_data = float(np.mean(y))
        print(
            f"noise abs={int(noise)} rms={int(noise_rms)} "
            f"5*mean={int(mean_data * 5)} 10*mean={int(mean_data * 10)}"
        )
    return noise
