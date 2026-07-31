import numpy as np
import pytest

from peer_n3_signature.core import (
    assert_background_equal,
    dominant_period,
    gaussian_covariance_ttee,
    weighted_project,
)


def test_weighted_project_removes_exact_derivative_span():
    x = np.linspace(-1.0, 1.0, 60)
    basis = np.column_stack([np.ones_like(x), x, x**2])
    signal = 2.0 * basis[:, 0] - 0.7 * basis[:, 1] + 0.2 * basis[:, 2]
    weights = np.linspace(1.0, 3.0, x.size)
    projected, coefficients = weighted_project(signal, basis, weights)
    assert np.linalg.norm(projected) < 1e-10
    assert np.allclose(coefficients, [2.0, -0.7, 0.2], atol=1e-10)


def test_gaussian_covariance_ttee_is_positive_definite():
    covariance = gaussian_covariance_ttee(
        ell=500,
        ctt=5.0e-4,
        cee=1.5e-5,
        cte=4.0e-5,
        fsky=0.6,
    )
    assert covariance.shape == (3, 3)
    assert np.all(np.linalg.eigvalsh(covariance) > 0)


def test_dominant_period_recovers_synthetic_oscillation():
    ell = np.arange(100, 1800)
    period_true = 275.0
    signal = np.sin(2 * np.pi * ell / period_true)
    recovered = dominant_period(ell, signal, min_period=100.0, max_period=500.0)
    assert recovered == pytest.approx(period_true, rel=0.03)


def test_background_gate_rejects_material_difference():
    reference = {
        "z": np.array([0.0, 10.0, 1000.0]),
        "H": np.array([70.0, 1400.0, 1.0e6]),
        "omega_de": np.array([0.7, 0.01, 0.08]),
        "derived": {"H0": 70.0, "rdrag": 143.0, "thetastar": 1.041},
    }
    trial = {
        "z": reference["z"].copy(),
        "H": reference["H"].copy(),
        "omega_de": reference["omega_de"].copy(),
        "derived": dict(reference["derived"]),
    }
    trial["H"][1] *= 1.001
    with pytest.raises(AssertionError, match="background mismatch"):
        assert_background_equal(reference, trial, rtol=1e-5, atol=1e-10)
