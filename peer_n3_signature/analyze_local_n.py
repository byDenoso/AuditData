from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .analyze import (
    cosine,
    derivative_models,
    feature_locations,
    fit_nuisance_coefficients,
    gaussian_covariance_ttee,
    load_model,
    normalized_channel,
    whiten_ttee,
)

CHANNEL_INDEX = {"TT": 0, "EE": 1, "TE": 3}
COV_INDEX = {"TT": 0, "EE": 1, "TE": 2}
BANDS = {
    "first_peaks": (30, 300),
    "planck_mid": (301, 800),
    "transition": (801, 1600),
    "damping": (1601, 2500),
    "deep_damping": (2501, 3000),
}
DELTA_N = 0.2


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


def channel_cv_norm(signal: np.ndarray, central: np.ndarray, channel: str, lo: int, hi: int) -> float:
    index = CHANNEL_INDEX[channel]
    covariance_index = COV_INDEX[channel]
    total = 0.0
    for ell in range(lo, hi + 1):
        covariance = gaussian_covariance_ttee(
            ell,
            float(central[ell, 0]),
            float(central[ell, 1]),
            float(central[ell, 3]),
            regularization=1e-12,
        )
        total += float(signal[ell, index] ** 2 / covariance[covariance_index, covariance_index])
    return float(np.sqrt(max(total, 0.0)))


def band_metrics(signal: np.ndarray, central: np.ndarray, lmax: int) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for name, (lo, hi) in BANDS.items():
        hi = min(hi, lmax)
        if hi < lo:
            continue
        combined = float(np.linalg.norm(whiten_ttee(signal, central, lo, hi)[0]))
        output[name] = {
            "ell_min": lo,
            "ell_max": hi,
            "combined_cv_norm": combined,
            "channels": {
                channel: channel_cv_norm(signal, central, channel, lo, hi)
                for channel in CHANNEL_INDEX
            },
        }
    return output


def fractional_summary(signal: np.ndarray, central: np.ndarray, lo: int, hi: int) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for channel in ("TT", "TE", "EE"):
        values = normalized_channel(signal, central, channel)[lo : hi + 1]
        values = values[np.isfinite(values)]
        output[channel] = {
            "rms_fraction": float(np.sqrt(np.mean(values**2))),
            "p95_abs_fraction": float(np.quantile(np.abs(values), 0.95)),
            "max_abs_fraction": float(np.max(np.abs(values))),
        }
    return output


def lens_cv_norm(signal: np.ndarray, central: np.ndarray, lo: int, hi: int) -> float:
    hi = min(hi, signal.shape[0] - 1, central.shape[0] - 1)
    multipoles = np.arange(lo, hi + 1)
    spectrum = central[multipoles, 0]
    variance = 2.0 * spectrum**2 / (2.0 * multipoles + 1.0)
    valid = np.isfinite(signal[multipoles, 0]) & np.isfinite(variance) & (variance > 0)
    return float(np.sqrt(np.sum(signal[multipoles[valid], 0] ** 2 / variance[valid])))


def write_csv(
    path: Path,
    ell: np.ndarray,
    local_step: np.ndarray,
    projected_step: np.ndarray,
    curvature: np.ndarray,
    projected_curvature: np.ndarray,
    central: np.ndarray,
    lens_step: np.ndarray,
    lens_projected: np.ndarray,
    central_lens: np.ndarray,
) -> None:
    normalized = {
        "local": {channel: normalized_channel(local_step, central, channel) for channel in CHANNEL_INDEX},
        "local_projected": {
            channel: normalized_channel(projected_step, central, channel) for channel in CHANNEL_INDEX
        },
        "curvature": {channel: normalized_channel(curvature, central, channel) for channel in CHANNEL_INDEX},
        "curvature_projected": {
            channel: normalized_channel(projected_curvature, central, channel) for channel in CHANNEL_INDEX
        },
    }
    phi_denominator = np.maximum(np.abs(central_lens[:, 0]), np.finfo(float).tiny)
    fields = ["ell"]
    for key in ("local", "local_projected", "curvature", "curvature_projected"):
        fields.extend(f"{key}_{channel}" for channel in ("TT", "TE", "EE"))
    fields.extend(("lens_local_step", "lens_local_step_projected"))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, multipole in enumerate(ell):
            row: dict[str, float | int] = {"ell": int(multipole)}
            for key in ("local", "local_projected", "curvature", "curvature_projected"):
                for channel in ("TT", "TE", "EE"):
                    row[f"{key}_{channel}"] = float(normalized[key][channel][index])
            row["lens_local_step"] = float(lens_step[index, 0] / phi_denominator[index])
            row["lens_local_step_projected"] = float(lens_projected[index, 0] / phi_denominator[index])
            writer.writerow(row)


def make_plots(
    output: Path,
    ell: np.ndarray,
    projected_step: np.ndarray,
    projected_curvature: np.ndarray,
    central: np.ndarray,
    bands: dict[str, Any],
    lens_step: np.ndarray,
    lens_projected: np.ndarray,
    central_lens: np.ndarray,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for axis, channel in zip(axes, ("TT", "TE", "EE")):
        axis.plot(ell[30:], 100.0 * normalized_channel(projected_step, central, channel)[30:])
        axis.axhline(0.0, linewidth=0.8)
        axis.set_ylabel(f"{channel} [%]")
    axes[-1].set_xlabel("multipole ell")
    axes[-1].set_xlim(30, int(ell[-1]))
    fig.suptitle("Projected local n=3 scalar fingerprint for Delta n=0.2")
    fig.tight_layout()
    fig.savefig(output / "local_n3_projected_signature.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for axis, channel in zip(axes, ("TT", "TE", "EE")):
        axis.plot(ell[30:], 100.0 * normalized_channel(projected_curvature, central, channel)[30:])
        axis.axhline(0.0, linewidth=0.8)
        axis.set_ylabel(f"{channel} [%]")
    axes[-1].set_xlabel("multipole ell")
    axes[-1].set_xlim(30, int(ell[-1]))
    fig.suptitle("n=3 deviation from linear interpolation between n=2.8 and n=3.2")
    fig.tight_layout()
    fig.savefig(output / "n3_curvature_signature.png", dpi=180)
    plt.close(fig)

    labels = list(bands)
    values = [bands[label]["combined_cv_norm"] for label in labels]
    fig, axis = plt.subplots(figsize=(11, 5))
    axis.bar(labels, values)
    axis.set_ylabel("cosmic-variance-limited norm")
    axis.set_title("Band ranking for the projected Delta n=0.2 template")
    axis.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(output / "local_n3_band_ranking.png", dpi=180)
    plt.close(fig)

    phi_denominator = np.maximum(np.abs(central_lens[:, 0]), np.finfo(float).tiny)
    fig, axis = plt.subplots(figsize=(12, 5))
    axis.plot(ell[2:], 100.0 * lens_step[2:, 0] / phi_denominator[2:], label="raw local n step")
    axis.plot(ell[2:], 100.0 * lens_projected[2:, 0] / phi_denominator[2:], label="after CMB nuisance projection")
    axis.axhline(0.0, linewidth=0.8)
    axis.set_xlim(2, min(2000, int(ell[-1])))
    axis.set_xlabel("lensing multipole L")
    axis.set_ylabel("Delta C_L^phiphi / C_L^phiphi [%]")
    axis.set_title("Lensing-potential companion template")
    axis.legend()
    fig.tight_layout()
    fig.savefig(output / "local_n3_lensing_signature.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    n28 = load_model(args.input, "scalar_n28")
    n30 = load_model(args.input, "scalar_n30")
    n32 = load_model(args.input, "scalar_n32")
    nuisance_names, derivatives = derivative_models(args.input)
    central = n30["lensed"]
    ell = np.asarray(n30["ell"], dtype=int)
    lmax = int(ell[-1])
    fit_max = min(3000, lmax)
    nuisance_cmb = [derivatives[name]["lensed"] for name in nuisance_names]

    dcdn = (n32["lensed"] - n28["lensed"]) / (2.0 * DELTA_N)
    local_step = DELTA_N * dcdn
    local_projected, coefficients, local_raw_norm, local_projected_norm = fit_nuisance_coefficients(
        local_step,
        nuisance_cmb,
        central,
        ell_min=30,
        ell_max=fit_max,
    )

    curvature = n30["lensed"] - 0.5 * (n28["lensed"] + n32["lensed"])
    curvature_projected, curvature_coefficients, curvature_raw_norm, curvature_projected_norm = (
        fit_nuisance_coefficients(
            curvature,
            nuisance_cmb,
            central,
            ell_min=30,
            ell_max=fit_max,
        )
    )

    local_white = whiten_ttee(local_projected, central, 30, fit_max)[0]
    curvature_white = whiten_ttee(curvature_projected, central, 30, fit_max)[0]
    curvature_cosine = cosine(local_white, curvature_white)

    lens_dcdn = (n32["lens"] - n28["lens"]) / (2.0 * DELTA_N)
    lens_step = DELTA_N * lens_dcdn
    lens_projected = lens_step - sum(
        coefficient * derivatives[name]["lens"]
        for coefficient, name in zip(coefficients, nuisance_names)
    )
    lensing = {
        "raw_cv_norm_L30_2000": lens_cv_norm(lens_step, n30["lens"], 30, 2000),
        "projected_cv_norm_L30_2000": lens_cv_norm(lens_projected, n30["lens"], 30, 2000),
        "projected_cv_norm_L30_400": lens_cv_norm(lens_projected, n30["lens"], 30, 400),
        "projected_cv_norm_L401_1000": lens_cv_norm(lens_projected, n30["lens"], 401, 1000),
        "projected_cv_norm_L1001_2000": lens_cv_norm(lens_projected, n30["lens"], 1001, 2000),
    }

    bands = band_metrics(local_projected, central, lmax)
    fractions = {
        name: fractional_summary(local_projected, central, values["ell_min"], values["ell_max"])
        for name, values in bands.items()
    }
    features = {
        channel: feature_locations(
            ell,
            normalized_channel(local_projected, central, channel),
            lo=100,
            hi=min(2500, lmax),
        )
        for channel in ("TT", "TE", "EE")
    }
    best_band_name, best_band = max(bands.items(), key=lambda item: item[1]["combined_cv_norm"])
    survival = local_projected_norm / local_raw_norm if local_raw_norm > 0 else 0.0
    curvature_ratio = curvature_projected_norm / local_projected_norm if local_projected_norm > 0 else 0.0

    if curvature_projected_norm >= 3.0:
        status = "N3_SPECIFIC_CURVATURE_PRESENT"
    elif local_projected_norm >= 2.0:
        status = "LOCAL_N3_SCALAR_FINGERPRINT_PRESENT_NOT_UNIQUE"
    else:
        status = "NO_RESOLVABLE_LOCAL_N3_FINGERPRINT"

    results = {
        "test_id": "T-N3-004",
        "status": status,
        "evidence_class": "PHYSICAL_LOCAL_N_DERIVATIVE_CAMB_FINGERPRINT",
        "model_family": "CAMB 1.6.6 EarlyQuintessence n=2.8, 3.0, 3.2",
        "central_n": 3.0,
        "delta_n_template": DELTA_N,
        "nuisance_parameters": nuisance_names,
        "nuisance_coefficients_local_step": {
            name: float(value) for name, value in zip(nuisance_names, coefficients)
        },
        "nuisance_coefficients_curvature": {
            name: float(value) for name, value in zip(nuisance_names, curvature_coefficients)
        },
        "local_step": {
            "raw_cv_norm_ell30_3000": local_raw_norm,
            "projected_cv_norm_ell30_3000": local_projected_norm,
            "survival_fraction": survival,
            "band_ranking": bands,
            "best_band": {"name": best_band_name, **best_band},
            "fractional_amplitudes": fractions,
            "features": features,
        },
        "n3_specific_curvature": {
            "definition": "C(n=3) - [C(n=2.8)+C(n=3.2)]/2",
            "raw_cv_norm_ell30_3000": curvature_raw_norm,
            "projected_cv_norm_ell30_3000": curvature_projected_norm,
            "ratio_to_local_step": curvature_ratio,
            "cosine_with_local_step": curvature_cosine,
        },
        "lensing_companion": lensing,
        "retired_test": {
            "test_id": "T-N3-003",
            "reason": (
                "Setting num_perturb_equations=0 preserves the background but produces a perturbatively "
                "non-closed counterfactual with pathological CMB amplitudes; it cannot identify a physical signature."
            ),
        },
        "claim_boundary": (
            "This is a physical local derivative within the scalar family and a cosmic-variance-limited template. "
            "It is not a Planck, ACT or SPT detection significance. A unique n=3 claim requires the curvature template, "
            "not the local derivative alone."
        ),
    }
    (args.output / "local_n3_metrics.json").write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_csv(
        args.output / "local_n3_spectra.csv",
        ell,
        local_step,
        local_projected,
        curvature,
        curvature_projected,
        central,
        lens_step,
        lens_projected,
        n30["lens"],
    )
    make_plots(
        args.output,
        ell,
        local_projected,
        curvature_projected,
        central,
        bands,
        lens_step,
        lens_projected,
        n30["lens"],
    )

    report = f"""# Physical local n=3 scalar CMB fingerprint

## Decision

**Status:** `{status}`

For a physical displacement **Delta n=0.2** around n=3, the nuisance-projected TT/TE/EE template has a cosmic-variance-limited norm of **{local_projected_norm:.3f}** over ell=30-{fit_max}. The raw norm is **{local_raw_norm:.3f}**, so **{100.0 * survival:.1f}%** survives projection of lnAs, ns, ombh2, omch2, tau, H0 and Alens.

The strongest band is **{best_band_name}**, ell={best_band['ell_min']}-{best_band['ell_max']}, with combined norm **{best_band['combined_cv_norm']:.3f}**. Its channel norms are TT={best_band['channels']['TT']:.3f}, TE={best_band['channels']['TE']:.3f}, and EE={best_band['channels']['EE']:.3f}.

## Morphology

- TT dominant period: {features['TT']['dominant_period']}
- TE dominant period: {features['TE']['dominant_period']}
- EE dominant period: {features['EE']['dominant_period']}

The n=3-specific curvature template has projected norm **{curvature_projected_norm:.3f}**, only **{100.0 * curvature_ratio:.1f}%** of the local-step norm. Therefore the measurable pattern is a scalar-family n-direction around n=3, not a unique cusp or isolated resonance at exactly n=3.

## Lensing companion

The projected lensing-potential ceiling is {lensing['projected_cv_norm_L30_2000']:.3f} over L=30-2000. Only {lensing['projected_cv_norm_L30_400']:.3f} lies at L<=400; most of the theoretical information is at L>1000.

## Scientific boundary

These norms assume ideal cosmic-variance-limited spectra. They identify where to build a fixed-template likelihood, not an observed significance. The next decisive test is a one-amplitude matched template fit in ACT DR6 and SPT-3G TE/EE, with Planck TT/TE/EE as the low-to-mid-ell anchor.
"""
    (args.output / "REPORT.md").write_text(report, encoding="utf-8")

    hashes = []
    for path in sorted(args.output.iterdir()):
        if path.is_file() and path.name != "files.sha256":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (args.output / "files.sha256").write_text("\n".join(hashes) + "\n", encoding="utf-8")

    print(
        "LOCAL_N3_SUMMARY="
        + json.dumps(
            json_safe(
                {
                    "status": status,
                    "projected_cv_norm": local_projected_norm,
                    "survival_fraction": survival,
                    "best_band": best_band_name,
                    "best_band_metrics": best_band,
                    "features": features,
                    "curvature_projected_cv_norm": curvature_projected_norm,
                    "lensing": lensing,
                }
            ),
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
