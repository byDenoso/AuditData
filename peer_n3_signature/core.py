from __future__ import annotations

from typing import Mapping

import numpy as np


def gaussian_covariance_ttee(
    ell: int,
    ctt: float,
    cee: float,
    cte: float,
    fsky: float = 1.0,
    regularization: float = 1e-14,
) -> np.ndarray:
    """Full-sky Gaussian covariance for (TT, EE, TE) at one multipole."""
    if ell < 2:
        raise ValueError("ell must be >= 2")
    if not 0 < fsky <= 1:
        raise ValueError("fsky must be in (0, 1]")
    nu = (2 * ell + 1) * fsky
    covariance = np.array(
        [
            [2 * ctt**2, 2 * cte**2, 2 * ctt * cte],
            [2 * cte**2, 2 * cee**2, 2 * cee * cte],
            [2 * ctt * cte, 2 * cee * cte, cte**2 + ctt * cee],
        ],
        dtype=float,
    ) / nu
    scale = max(float(np.max(np.abs(np.diag(covariance)))), np.finfo(float).tiny)
    covariance += np.eye(3) * regularization * scale
    return covariance


def weighted_project(
    signal: np.ndarray,
    basis: np.ndarray,
    weight_or_precision: np.ndarray,
    rcond: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """Remove the weighted least-squares span of basis columns from signal."""
    y = np.asarray(signal, dtype=float)
    design = np.asarray(basis, dtype=float)
    weight = np.asarray(weight_or_precision, dtype=float)
    if y.ndim != 1 or design.ndim != 2 or design.shape[0] != y.size:
        raise ValueError("incompatible signal and basis shapes")
    if weight.ndim == 1:
        if weight.size != y.size or np.any(weight <= 0):
            raise ValueError("weights must be positive and match signal length")
        sqrt_w = np.sqrt(weight)
        lhs = design * sqrt_w[:, None]
        rhs = y * sqrt_w
    elif weight.ndim == 2:
        if weight.shape != (y.size, y.size):
            raise ValueError("precision matrix has wrong shape")
        eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (weight + weight.T))
        floor = max(np.max(eigenvalues) * rcond, np.finfo(float).tiny)
        sqrt_precision = (eigenvectors * np.sqrt(np.maximum(eigenvalues, floor))) @ eigenvectors.T
        lhs = sqrt_precision @ design
        rhs = sqrt_precision @ y
    else:
        raise ValueError("weight_or_precision must be 1D or 2D")
    coefficients, *_ = np.linalg.lstsq(lhs, rhs, rcond=rcond)
    residual = y - design @ coefficients
    return residual, coefficients


def dominant_period(
    ell: np.ndarray,
    signal: np.ndarray,
    min_period: float,
    max_period: float,
) -> float:
    """Estimate dominant oscillation period on a regularly sampled ell grid."""
    x = np.asarray(ell, dtype=float)
    y = np.asarray(signal, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size or x.size < 8:
        raise ValueError("ell and signal must be matching 1D arrays")
    spacing = np.diff(x)
    if not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-10):
        raise ValueError("ell grid must be regular")
    if not 0 < min_period < max_period:
        raise ValueError("invalid period bounds")
    trend = np.polyval(np.polyfit(x, y, deg=min(2, x.size - 1)), x)
    centered = y - trend
    window = np.hanning(x.size)
    spectrum = np.abs(np.fft.rfft(centered * window)) ** 2
    frequencies = np.fft.rfftfreq(x.size, d=spacing[0])
    periods = np.full_like(frequencies, np.inf)
    periods[1:] = 1.0 / frequencies[1:]
    valid = (periods >= min_period) & (periods <= max_period)
    if not np.any(valid):
        raise ValueError("period bounds exclude all Fourier modes")
    valid_indices = np.flatnonzero(valid)
    peak = valid_indices[int(np.argmax(spectrum[valid]))]
    refined_bin = float(peak)
    if 0 < peak < spectrum.size - 1:
        left, center, right = np.log(np.maximum(spectrum[peak - 1 : peak + 2], np.finfo(float).tiny))
        curvature = left - 2.0 * center + right
        if abs(curvature) > 1e-15:
            refined_bin += 0.5 * (left - right) / curvature
    refined_frequency = refined_bin / (x.size * spacing[0])
    return float(1.0 / refined_frequency)


def assert_background_equal(
    reference: Mapping[str, object],
    trial: Mapping[str, object],
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> dict[str, float]:
    """Fail closed when exact-background perturbation controls drift materially."""
    metrics: dict[str, float] = {}
    for key in ("z", "H", "omega_de"):
        ref = np.asarray(reference[key], dtype=float)
        got = np.asarray(trial[key], dtype=float)
        if ref.shape != got.shape:
            raise AssertionError(f"background mismatch: {key} shape {ref.shape} != {got.shape}")
        denom = np.maximum(np.abs(ref), atol)
        metrics[f"max_rel_{key}"] = float(np.max(np.abs(got - ref) / denom))
        if not np.allclose(ref, got, rtol=rtol, atol=atol):
            raise AssertionError(f"background mismatch: {key} max_rel={metrics[f'max_rel_{key}']:.6g}")
    ref_derived = reference["derived"]
    trial_derived = trial["derived"]
    for key in ("H0", "rdrag", "thetastar"):
        ref = float(ref_derived[key])
        got = float(trial_derived[key])
        rel = abs(got - ref) / max(abs(ref), atol)
        metrics[f"rel_{key}"] = rel
        if not np.isclose(ref, got, rtol=rtol, atol=atol):
            raise AssertionError(f"background mismatch: {key} rel={rel:.6g}")
    return metrics
