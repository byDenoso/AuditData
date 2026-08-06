from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from .core import assert_background_equal, dominant_period, gaussian_covariance_ttee
from .generate import DERIVATIVE_STEPS

CHANNEL_COLUMNS = {"TT": 0, "EE": 1, "TE": 3}
BANDS = {
    "low": (2, 29),
    "first_peaks": (30, 300),
    "planck_mid": (301, 800),
    "transition": (801, 1600),
    "damping": (1601, 2500),
    "deep_damping": (2501, 3000),
}


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def load_model(root: Path, label: str) -> dict[str, Any]:
    arrays = np.load(root / f"{label}.npz")
    metadata = json.loads((root / f"{label}.json").read_text(encoding="utf-8"))
    return {key: np.asarray(arrays[key]) for key in arrays.files} | {"metadata": metadata}


def background_payload(model: dict[str, Any]) -> dict[str, Any]:
    return {
        "z": model["z"],
        "H": model["H"],
        "omega_de": model["omega_de"],
        "derived": model["metadata"]["derived"],
    }


def derivative_models(root: Path) -> tuple[list[str], dict[str, dict[str, np.ndarray]]]:
    names = list(DERIVATIVE_STEPS)
    output: dict[str, dict[str, np.ndarray]] = {}
    for name in names:
        minus = load_model(root, f"deriv_{name}_minus")
        plus = load_model(root, f"deriv_{name}_plus")
        scale = 2.0 * DERIVATIVE_STEPS[name]
        output[name] = {
            "lensed": (plus["lensed"] - minus["lensed"]) / scale,
            "unlensed": (plus["unlensed"] - minus["unlensed"]) / scale,
            "lens": (plus["lens"] - minus["lens"]) / scale,
        }
    return names, output


def whiten_ttee(
    signal: np.ndarray,
    central: np.ndarray,
    ell_min: int,
    ell_max: int,
    fsky: float = 1.0,
) -> tuple[np.ndarray, list[np.ndarray]]:
    rows: list[np.ndarray] = []
    cholesky: list[np.ndarray] = []
    for ell in range(ell_min, ell_max + 1):
        covariance = gaussian_covariance_ttee(
            ell,
            ctt=float(central[ell, 0]),
            cee=float(central[ell, 1]),
            cte=float(central[ell, 3]),
            fsky=fsky,
            regularization=1e-12,
        )
        root = np.linalg.cholesky(covariance)
        vector = signal[ell, [0, 1, 3]]
        rows.append(np.linalg.solve(root, vector))
        cholesky.append(root)
    return np.concatenate(rows), cholesky


def whiten_ttee_basis(
    derivatives: list[np.ndarray],
    central: np.ndarray,
    ell_min: int,
    ell_max: int,
    fsky: float = 1.0,
) -> np.ndarray:
    columns = []
    for derivative in derivatives:
        whitened, _ = whiten_ttee(derivative, central, ell_min, ell_max, fsky=fsky)
        columns.append(whitened)
    return np.column_stack(columns)


def fit_nuisance_coefficients(
    signal: np.ndarray,
    derivative_arrays: list[np.ndarray],
    central: np.ndarray,
    ell_min: int,
    ell_max: int,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    whitened_signal, _ = whiten_ttee(signal, central, ell_min, ell_max)
    whitened_basis = whiten_ttee_basis(derivative_arrays, central, ell_min, ell_max)
    coefficients, *_ = np.linalg.lstsq(whitened_basis, whitened_signal, rcond=1e-10)
    projected = signal - sum(coefficient * derivative for coefficient, derivative in zip(coefficients, derivative_arrays))
    whitened_projected, _ = whiten_ttee(projected, central, ell_min, ell_max)
    return projected, coefficients, float(np.linalg.norm(whitened_signal)), float(np.linalg.norm(whitened_projected))


def band_norms(
    raw: np.ndarray,
    projected: np.ndarray,
    central: np.ndarray,
    lmax: int,
) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for band, (lo, hi) in BANDS.items():
        hi = min(hi, lmax)
        if hi < lo:
            continue
        raw_white, _ = whiten_ttee(raw, central, lo, hi)
        projected_white, _ = whiten_ttee(projected, central, lo, hi)
        channel_metrics: dict[str, dict[str, float]] = {}
        for channel, index in CHANNEL_COLUMNS.items():
            raw_snr2 = 0.0
            projected_snr2 = 0.0
            covariance_index = {"TT": 0, "EE": 1, "TE": 2}[channel]
            for ell in range(lo, hi + 1):
                covariance = gaussian_covariance_ttee(
                    ell,
                    float(central[ell, 0]),
                    float(central[ell, 1]),
                    float(central[ell, 3]),
                    regularization=1e-12,
                )
                variance = covariance[covariance_index, covariance_index]
                raw_snr2 += float(raw[ell, index] ** 2 / variance)
                projected_snr2 += float(projected[ell, index] ** 2 / variance)
            channel_metrics[channel] = {
                "raw_cv_snr": float(np.sqrt(max(raw_snr2, 0.0))),
                "projected_cv_snr": float(np.sqrt(max(projected_snr2, 0.0))),
            }
        output[band] = {
            "ell_min": lo,
            "ell_max": hi,
            "combined_raw_cv_snr": float(np.linalg.norm(raw_white)),
            "combined_projected_cv_snr": float(np.linalg.norm(projected_white)),
            "channels": channel_metrics,
        }
    return output


def lensing_projection(
    signal: np.ndarray,
    derivative_arrays: list[np.ndarray],
    coefficients: np.ndarray,
    central: np.ndarray,
    ell_min: int,
    ell_max: int,
) -> dict[str, Any]:
    projected = signal - sum(coefficient * derivative for coefficient, derivative in zip(coefficients, derivative_arrays))
    central_phi = central[:, 0]
    valid = np.arange(ell_min, ell_max + 1)
    variance = 2.0 * np.maximum(central_phi[valid] ** 2, np.finfo(float).tiny) / (2.0 * valid + 1.0)
    raw_snr = float(np.sqrt(np.sum(signal[valid, 0] ** 2 / variance)))
    projected_snr = float(np.sqrt(np.sum(projected[valid, 0] ** 2 / variance)))
    return {"projected": projected, "raw_cv_snr": raw_snr, "projected_cv_snr": projected_snr}


def normalized_channel(signal: np.ndarray, central: np.ndarray, channel: str) -> np.ndarray:
    if channel == "TT":
        denominator = np.maximum(np.abs(central[:, 0]), np.finfo(float).tiny)
        return signal[:, 0] / denominator
    if channel == "EE":
        denominator = np.maximum(np.abs(central[:, 1]), np.finfo(float).tiny)
        return signal[:, 1] / denominator
    if channel == "TE":
        denominator = np.sqrt(np.maximum(np.abs(central[:, 0] * central[:, 1]), np.finfo(float).tiny))
        return signal[:, 3] / denominator
    raise ValueError(channel)


def feature_locations(ell: np.ndarray, normalized: np.ndarray, lo: int = 100, hi: int = 2500) -> dict[str, Any]:
    selection = (ell >= lo) & (ell <= hi) & np.isfinite(normalized)
    x = ell[selection]
    y = normalized[selection]
    if x.size < 20 or np.max(np.abs(y)) == 0:
        return {"dominant_period": None, "zero_crossings": [], "peaks": []}
    period = dominant_period(x, y, min_period=80.0, max_period=600.0)
    signs = np.signbit(y)
    zeros = x[1:][signs[1:] != signs[:-1]]
    peak_indices, _ = find_peaks(np.abs(y), distance=40, prominence=max(np.max(np.abs(y)) * 0.03, 1e-14))
    ranking = sorted(peak_indices, key=lambda index: abs(y[index]), reverse=True)[:12]
    peaks = [{"ell": int(x[index]), "amplitude": float(y[index])} for index in ranking]
    return {
        "dominant_period": float(period),
        "zero_crossings": [int(value) for value in zeros[:30]],
        "peaks": peaks,
    }


def cosine(first: np.ndarray, second: np.ndarray) -> float:
    denominator = np.linalg.norm(first) * np.linalg.norm(second)
    return float(np.dot(first, second) / denominator) if denominator > 0 else float("nan")


def write_csv(
    path: Path,
    ell: np.ndarray,
    raw: np.ndarray,
    projected: np.ndarray,
    fluid: np.ndarray,
    central: np.ndarray,
    lens_raw: np.ndarray,
    lens_projected: np.ndarray,
    central_lens: np.ndarray,
) -> None:
    fields = [
        "ell",
        "raw_TT_over_TT",
        "raw_EE_over_EE",
        "raw_TE_over_sqrt_TT_EE",
        "projected_TT_over_TT",
        "projected_EE_over_EE",
        "projected_TE_over_sqrt_TT_EE",
        "scalar_minus_fluid_TT_over_TT",
        "scalar_minus_fluid_EE_over_EE",
        "scalar_minus_fluid_TE_over_sqrt_TT_EE",
        "raw_phi_over_phi",
        "projected_phi_over_phi",
    ]
    normalizations = {
        "raw": {channel: normalized_channel(raw, central, channel) for channel in CHANNEL_COLUMNS},
        "projected": {channel: normalized_channel(projected, central, channel) for channel in CHANNEL_COLUMNS},
        "fluid": {channel: normalized_channel(fluid, central, channel) for channel in CHANNEL_COLUMNS},
    }
    phi_denominator = np.maximum(np.abs(central_lens[:, 0]), np.finfo(float).tiny)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, value in enumerate(ell):
            writer.writerow(
                {
                    "ell": int(value),
                    "raw_TT_over_TT": normalizations["raw"]["TT"][index],
                    "raw_EE_over_EE": normalizations["raw"]["EE"][index],
                    "raw_TE_over_sqrt_TT_EE": normalizations["raw"]["TE"][index],
                    "projected_TT_over_TT": normalizations["projected"]["TT"][index],
                    "projected_EE_over_EE": normalizations["projected"]["EE"][index],
                    "projected_TE_over_sqrt_TT_EE": normalizations["projected"]["TE"][index],
                    "scalar_minus_fluid_TT_over_TT": normalizations["fluid"]["TT"][index],
                    "scalar_minus_fluid_EE_over_EE": normalizations["fluid"]["EE"][index],
                    "scalar_minus_fluid_TE_over_sqrt_TT_EE": normalizations["fluid"]["TE"][index],
                    "raw_phi_over_phi": lens_raw[index, 0] / phi_denominator[index],
                    "projected_phi_over_phi": lens_projected[index, 0] / phi_denominator[index],
                }
            )


def make_plots(
    output: Path,
    ell: np.ndarray,
    raw: np.ndarray,
    projected: np.ndarray,
    fluid: np.ndarray,
    central: np.ndarray,
    bands: dict[str, dict[str, Any]],
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for axis, channel in zip(axes, ("TT", "TE", "EE")):
        axis.plot(ell[2:], normalized_channel(raw, central, channel)[2:], label="scalar n3 - same background, no perturbations")
        axis.axhline(0.0, linewidth=0.8)
        axis.set_ylabel(channel)
        axis.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("multipole ell")
    axes[-1].set_xlim(2, int(ell[-1]))
    fig.suptitle("Exact-background scalar perturbation signature")
    fig.tight_layout()
    fig.savefig(output / "exact_scalar_signature.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for axis, channel in zip(axes, ("TT", "TE", "EE")):
        axis.plot(ell[2:], normalized_channel(projected, central, channel)[2:], label="after nuisance projection")
        axis.plot(ell[2:], normalized_channel(fluid, central, channel)[2:], alpha=0.7, label="scalar n3 - matched effective fluid")
        axis.axhline(0.0, linewidth=0.8)
        axis.set_ylabel(channel)
        axis.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("multipole ell")
    axes[-1].set_xlim(2, int(ell[-1]))
    fig.suptitle("Surviving n3 fingerprint and fluid-control residual")
    fig.tight_layout()
    fig.savefig(output / "projected_signature.png", dpi=170)
    plt.close(fig)

    labels = list(bands)
    values = [bands[label]["combined_projected_cv_snr"] for label in labels]
    fig, axis = plt.subplots(figsize=(11, 5))
    axis.bar(labels, values)
    axis.set_ylabel("cosmic-variance-limited projected norm")
    axis.set_title("Multipole-band ranking of the surviving signature")
    axis.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(output / "band_ranking.png", dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    scalar = load_model(args.input, "scalar_n30")
    nopert = load_model(args.input, "scalar_n30_nopert")
    fluid = load_model(args.input, "fluid_matched")
    n28 = load_model(args.input, "scalar_n28")
    n32 = load_model(args.input, "scalar_n32")
    names, derivatives = derivative_models(args.input)

    background_metrics = assert_background_equal(
        background_payload(scalar),
        background_payload(nopert),
        rtol=2e-9,
        atol=1e-10,
    )
    ell = scalar["ell"]
    lmax = int(ell[-1])
    fit_max = min(3000, lmax)
    exact_raw = scalar["lensed"] - nopert["lensed"]
    exact_unlensed = scalar["unlensed"] - nopert["unlensed"]
    exact_lens = scalar["lens"] - nopert["lens"]
    fluid_residual = scalar["lensed"] - fluid["lensed"]
    derivative_lensed = [derivatives[name]["lensed"] for name in names]
    derivative_unlensed = [derivatives[name]["unlensed"] for name in names]
    derivative_lens = [derivatives[name]["lens"] for name in names]

    projected, coefficients, raw_norm, projected_norm = fit_nuisance_coefficients(
        exact_raw,
        derivative_lensed,
        scalar["lensed"],
        ell_min=30,
        ell_max=fit_max,
    )
    projected_unlensed = exact_unlensed - sum(
        coefficient * derivative for coefficient, derivative in zip(coefficients, derivative_unlensed)
    )
    lens = lensing_projection(
        exact_lens,
        derivative_lens,
        coefficients,
        scalar["lens"],
        ell_min=30,
        ell_max=fit_max,
    )
    bands = band_norms(exact_raw, projected, scalar["lensed"], lmax)

    exact_white, _ = whiten_ttee(projected, scalar["lensed"], 30, fit_max)
    n_derivative = (n32["lensed"] - n28["lensed"]) / 0.4
    n_projected, _, _, _ = fit_nuisance_coefficients(
        n_derivative,
        derivative_lensed,
        scalar["lensed"],
        ell_min=30,
        ell_max=fit_max,
    )
    n_white, _ = whiten_ttee(n_projected, scalar["lensed"], 30, fit_max)
    fluid_projected, _, _, _ = fit_nuisance_coefficients(
        fluid_residual,
        derivative_lensed,
        scalar["lensed"],
        ell_min=30,
        ell_max=fit_max,
    )
    fluid_white, _ = whiten_ttee(fluid_projected, scalar["lensed"], 30, fit_max)

    features = {
        channel: feature_locations(ell, normalized_channel(projected, scalar["lensed"], channel))
        for channel in CHANNEL_COLUMNS
    }
    best_band_name, best_band = max(bands.items(), key=lambda item: item[1]["combined_projected_cv_snr"])
    survival = projected_norm / raw_norm if raw_norm > 0 else 0.0
    n_cosine = cosine(exact_white, n_white)
    fluid_cosine = cosine(exact_white, fluid_white)
    if projected_norm >= 5.0 and survival >= 0.20 and abs(n_cosine) >= 0.15:
        status = "DISTINCT_N3_SCALAR_PERTURBATION_FINGERPRINT"
    elif projected_norm >= 3.0:
        status = "SCALAR_PERTURBATION_RESPONSE_NOT_N3_SPECIFIC"
    else:
        status = "NO_DISTINCT_N3_SCALAR_PERTURBATION_FINGERPRINT"

    fluid_match = json.loads((args.input / "fluid_match.json").read_text(encoding="utf-8"))
    results = {
        "test_id": "T-N3-003",
        "status": status,
        "evidence_class": "EXACT_BACKGROUND_CAMB_PERTURBATION_ABLATION",
        "configuration": json.loads((args.input / "generation_contract.json").read_text(encoding="utf-8")),
        "background_equality": background_metrics,
        "nuisance_parameters": names,
        "nuisance_coefficients": {name: float(value) for name, value in zip(names, coefficients)},
        "raw_cv_norm_ell30_3000": raw_norm,
        "projected_cv_norm_ell30_3000": projected_norm,
        "survival_fraction": survival,
        "lensing": {key: value for key, value in lens.items() if key != "projected"},
        "band_ranking": bands,
        "best_band": {"name": best_band_name, **best_band},
        "features": features,
        "cosine_with_local_n_derivative": n_cosine,
        "cosine_with_scalar_minus_matched_fluid": fluid_cosine,
        "fluid_background_match": fluid_match,
        "claim_boundary": (
            "This is a theoretical CAMB fingerprint from an exact-background perturbation ablation. "
            "Cosmic-variance norms are detectability ceilings, not Planck, ACT or SPT likelihood significances."
        ),
    }
    (args.output / "signature_metrics.json").write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True), encoding="utf-8"
    )
    write_csv(
        args.output / "signature_spectra.csv",
        ell,
        exact_raw,
        projected,
        fluid_residual,
        scalar["lensed"],
        exact_lens,
        lens["projected"],
        scalar["lens"],
    )
    make_plots(args.output, ell, exact_raw, projected, fluid_residual, scalar["lensed"], bands)

    report = f"""# PEER n=3 scalar perturbation signature gate

## Result

**Status:** `{status}`

The exact-background ablation leaves a projected cosmic-variance-limited norm of **{projected_norm:.3f}** over ell=30-{fit_max}, compared with **{raw_norm:.3f}** before projecting the standard cosmological directions. The surviving fraction is **{survival:.3f}**.

The strongest multipole band is **{best_band_name}** (ell={best_band['ell_min']}-{best_band['ell_max']}), with projected norm **{best_band['combined_projected_cv_snr']:.3f}**.

Correlation with the local n derivative is **{n_cosine:.3f}**. Correlation with the scalar-minus-effective-fluid residual is **{fluid_cosine:.3f}**.

## Background gate

The scalar-full and scalar-no-perturbation builds use the same EarlyQuintessence background. Maximum relative differences:

- H(z): {background_metrics['max_rel_H']:.3e}
- Omega_PEER(z): {background_metrics['max_rel_omega_de']:.3e}
- H0: {background_metrics['rel_H0']:.3e}
- rdrag: {background_metrics['rel_rdrag']:.3e}
- thetastar: {background_metrics['rel_thetastar']:.3e}

## Spectral morphology

- TT dominant period: {features['TT']['dominant_period']}
- TE dominant period: {features['TE']['dominant_period']}
- EE dominant period: {features['EE']['dominant_period']}

## Interpretation boundary

This result identifies or rejects a theoretical linear-perturbation fingerprint of the n=3 scalar under the canonical PEER slice. It is not an observed detection, a matched likelihood comparison, or a model-selection result. Actual Planck/ACT/SPT significance requires convolving this fixed template with each experiment's windows, covariance and nuisance model.
"""
    (args.output / "REPORT.md").write_text(report, encoding="utf-8")

    hashes = []
    for path in sorted(args.output.iterdir()):
        if path.is_file() and path.name != "files.sha256":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (args.output / "files.sha256").write_text("\n".join(hashes) + "\n", encoding="utf-8")
    print("N3_SIGNATURE_SUMMARY=" + json.dumps(json_safe({
        "status": status,
        "projected_cv_norm": projected_norm,
        "survival_fraction": survival,
        "best_band": best_band_name,
        "n_cosine": n_cosine,
        "fluid_cosine": fluid_cosine,
        "features": features,
    }), sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
