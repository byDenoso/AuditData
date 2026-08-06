from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .analyze import (
    band_norms,
    background_payload,
    cosine,
    derivative_models,
    feature_locations,
    fit_nuisance_coefficients,
    lensing_projection,
    load_model,
    make_plots,
    normalized_channel,
    whiten_ttee,
)
from .core import assert_background_equal

CHANNELS = ("TT", "TE", "EE")


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


def write_signature_csv(
    path: Path,
    ell: np.ndarray,
    exact_raw: np.ndarray,
    exact_projected: np.ndarray,
    fluid_raw: np.ndarray,
    fluid_projected: np.ndarray,
    central_cmb: np.ndarray,
    lens_raw: np.ndarray,
    lens_projected: np.ndarray,
    central_lens: np.ndarray,
) -> None:
    exact_norm = {channel: normalized_channel(exact_raw, central_cmb, channel) for channel in CHANNELS}
    exact_proj_norm = {channel: normalized_channel(exact_projected, central_cmb, channel) for channel in CHANNELS}
    fluid_norm = {channel: normalized_channel(fluid_raw, central_cmb, channel) for channel in CHANNELS}
    fluid_proj_norm = {channel: normalized_channel(fluid_projected, central_cmb, channel) for channel in CHANNELS}
    phi_denominator = np.maximum(np.abs(central_lens[:, 0]), np.finfo(float).tiny)
    fieldnames = [
        "ell",
        "exact_raw_TT",
        "exact_raw_TE",
        "exact_raw_EE",
        "exact_projected_TT",
        "exact_projected_TE",
        "exact_projected_EE",
        "fluid_raw_TT",
        "fluid_raw_TE",
        "fluid_raw_EE",
        "fluid_projected_TT",
        "fluid_projected_TE",
        "fluid_projected_EE",
        "fluid_lensing_raw_phi",
        "fluid_lensing_projected_phi",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, multipole in enumerate(ell):
            writer.writerow(
                {
                    "ell": int(multipole),
                    "exact_raw_TT": exact_norm["TT"][index],
                    "exact_raw_TE": exact_norm["TE"][index],
                    "exact_raw_EE": exact_norm["EE"][index],
                    "exact_projected_TT": exact_proj_norm["TT"][index],
                    "exact_projected_TE": exact_proj_norm["TE"][index],
                    "exact_projected_EE": exact_proj_norm["EE"][index],
                    "fluid_raw_TT": fluid_norm["TT"][index],
                    "fluid_raw_TE": fluid_norm["TE"][index],
                    "fluid_raw_EE": fluid_norm["EE"][index],
                    "fluid_projected_TT": fluid_proj_norm["TT"][index],
                    "fluid_projected_TE": fluid_proj_norm["TE"][index],
                    "fluid_projected_EE": fluid_proj_norm["EE"][index],
                    "fluid_lensing_raw_phi": lens_raw[index, 0] / phi_denominator[index],
                    "fluid_lensing_projected_phi": lens_projected[index, 0] / phi_denominator[index],
                }
            )


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
    nuisance_names, derivatives = derivative_models(args.input)

    background_metrics = assert_background_equal(
        background_payload(scalar),
        background_payload(nopert),
        rtol=2e-9,
        atol=1e-10,
    )
    if nopert["metadata"].get("do_lensing", True):
        raise RuntimeError("exact no-perturbation control must have lensing disabled")

    ell = np.asarray(scalar["ell"], dtype=int)
    lmax = int(ell[-1])
    fit_max = min(3000, lmax)
    central_cmb = scalar["unlensed"]
    derivative_cmb = [derivatives[name]["unlensed"] for name in nuisance_names]

    exact_raw = scalar["unlensed"] - nopert["unlensed"]
    exact_projected, exact_coefficients, exact_raw_norm, exact_projected_norm = fit_nuisance_coefficients(
        exact_raw,
        derivative_cmb,
        central_cmb,
        ell_min=30,
        ell_max=fit_max,
    )

    fluid_raw = scalar["unlensed"] - fluid["unlensed"]
    fluid_projected, fluid_coefficients, fluid_raw_norm, fluid_projected_norm = fit_nuisance_coefficients(
        fluid_raw,
        derivative_cmb,
        central_cmb,
        ell_min=30,
        ell_max=fit_max,
    )

    n_derivative = (n32["unlensed"] - n28["unlensed"]) / 0.4
    n_projected, _, _, _ = fit_nuisance_coefficients(
        n_derivative,
        derivative_cmb,
        central_cmb,
        ell_min=30,
        ell_max=fit_max,
    )

    exact_white, _ = whiten_ttee(exact_projected, central_cmb, 30, fit_max)
    fluid_white, _ = whiten_ttee(fluid_projected, central_cmb, 30, fit_max)
    n_white, _ = whiten_ttee(n_projected, central_cmb, 30, fit_max)
    n_cosine = cosine(exact_white, n_white)
    fluid_cosine = cosine(exact_white, fluid_white)

    lens_raw = scalar["lens"] - fluid["lens"]
    lens_result = lensing_projection(
        lens_raw,
        [derivatives[name]["lens"] for name in nuisance_names],
        fluid_coefficients,
        scalar["lens"],
        ell_min=30,
        ell_max=fit_max,
    )

    bands = band_norms(exact_raw, exact_projected, central_cmb, lmax)
    features = {
        channel: feature_locations(
            ell,
            normalized_channel(exact_projected, central_cmb, channel),
        )
        for channel in CHANNELS
    }
    best_band_name, best_band = max(
        bands.items(),
        key=lambda item: item[1]["combined_projected_cv_snr"],
    )
    survival = exact_projected_norm / exact_raw_norm if exact_raw_norm > 0 else 0.0

    if exact_projected_norm >= 5.0 and survival >= 0.20 and abs(n_cosine) >= 0.15:
        status = "DISTINCT_N3_SCALAR_PERTURBATION_FINGERPRINT"
    elif exact_projected_norm >= 3.0:
        status = "SCALAR_PERTURBATION_RESPONSE_NOT_N3_SPECIFIC"
    else:
        status = "NO_DISTINCT_N3_SCALAR_PERTURBATION_FINGERPRINT"

    fluid_match = json.loads((args.input / "fluid_match.json").read_text(encoding="utf-8"))
    results = {
        "test_id": "T-N3-003",
        "status": status,
        "evidence_class": "EXACT_BACKGROUND_CAMB_PERTURBATION_ABLATION",
        "configuration": json.loads((args.input / "generation_contract.json").read_text(encoding="utf-8")),
        "exact_ablation_space": "unlensed TT/TE/EE",
        "background_equality": background_metrics,
        "nuisance_parameters": nuisance_names,
        "nuisance_coefficients_exact": {
            name: float(value) for name, value in zip(nuisance_names, exact_coefficients)
        },
        "nuisance_coefficients_fluid_control": {
            name: float(value) for name, value in zip(nuisance_names, fluid_coefficients)
        },
        "raw_cv_norm_ell30_3000": exact_raw_norm,
        "projected_cv_norm_ell30_3000": exact_projected_norm,
        "survival_fraction": survival,
        "fluid_control": {
            "raw_cv_norm_ell30_3000": fluid_raw_norm,
            "projected_cv_norm_ell30_3000": fluid_projected_norm,
            "background_match": fluid_match,
        },
        "lensing_control": {
            "type": "scalar_n3_minus_background_matched_effective_fluid",
            "raw_cv_snr": lens_result["raw_cv_snr"],
            "projected_cv_snr": lens_result["projected_cv_snr"],
            "exact_no_perturbation_lensing_available": False,
        },
        "band_ranking": bands,
        "best_band": {"name": best_band_name, **best_band},
        "features": features,
        "cosine_with_local_n_derivative": n_cosine,
        "cosine_with_scalar_minus_matched_fluid": fluid_cosine,
        "claim_boundary": (
            "This is a theoretical CAMB fingerprint from an exact-background perturbation ablation in unlensed TT/TE/EE. "
            "Lensing is reported only as a matched-fluid control. Cosmic-variance norms are detectability ceilings, "
            "not Planck, ACT or SPT likelihood significances."
        ),
    }
    (args.output / "signature_metrics.json").write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True),
        encoding="utf-8",
    )

    write_signature_csv(
        args.output / "signature_spectra.csv",
        ell,
        exact_raw,
        exact_projected,
        fluid_raw,
        fluid_projected,
        central_cmb,
        lens_raw,
        lens_result["projected"],
        scalar["lens"],
    )
    make_plots(
        args.output,
        ell,
        exact_raw,
        exact_projected,
        fluid_raw,
        central_cmb,
        bands,
    )

    report = f"""# PEER n=3 scalar perturbation signature gate

## Result

**Status:** `{status}`

The exact-background ablation in **unlensed TT/TE/EE** leaves a projected cosmic-variance-limited norm of **{exact_projected_norm:.3f}** over ell=30-{fit_max}, compared with **{exact_raw_norm:.3f}** before projecting the standard cosmological directions. The surviving fraction is **{survival:.3f}**.

The strongest multipole band is **{best_band_name}** (ell={best_band['ell_min']}-{best_band['ell_max']}), with projected norm **{best_band['combined_projected_cv_snr']:.3f}**.

Correlation with the local n derivative is **{n_cosine:.3f}**. Correlation with the scalar-minus-effective-fluid residual is **{fluid_cosine:.3f}**.

## Background gate

The scalar-full and scalar-no-perturbation builds share the same EarlyQuintessence background. Maximum relative differences:

- H(z): {background_metrics['max_rel_H']:.3e}
- Omega_PEER(z): {background_metrics['max_rel_omega_de']:.3e}
- H0: {background_metrics['rel_H0']:.3e}
- rdrag: {background_metrics['rel_rdrag']:.3e}
- thetastar: {background_metrics['rel_thetastar']:.3e}

## Spectral morphology

- TT dominant period: {features['TT']['dominant_period']}
- TE dominant period: {features['TE']['dominant_period']}
- EE dominant period: {features['EE']['dominant_period']}

## Lensing control

An exact no-perturbation lensed counterfactual is not defined by CAMB because removing the scalar perturbation equations removes a self-consistent lensing normalization. The lensing result therefore compares the scalar n=3 model with the background-matched effective-fluid control. Its projected cosmic-variance-limited norm is **{lens_result['projected_cv_snr']:.3f}**.

## Interpretation boundary

This result identifies or rejects a theoretical linear-perturbation fingerprint of the n=3 scalar under the canonical PEER slice. It is not an observed detection, a matched likelihood comparison, or a model-selection result. Actual Planck/ACT/SPT significance requires convolving this fixed template with each experiment's windows, covariance and nuisance model.
"""
    (args.output / "REPORT.md").write_text(report, encoding="utf-8")

    hashes = []
    for path in sorted(args.output.iterdir()):
        if path.is_file() and path.name != "files.sha256":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (args.output / "files.sha256").write_text("\n".join(hashes) + "\n", encoding="utf-8")

    summary = {
        "status": status,
        "projected_cv_norm": exact_projected_norm,
        "survival_fraction": survival,
        "best_band": best_band_name,
        "n_cosine": n_cosine,
        "fluid_cosine": fluid_cosine,
        "features": features,
        "lensing_control_projected_cv_snr": lens_result["projected_cv_snr"],
    }
    print("N3_SIGNATURE_SUMMARY=" + json.dumps(json_safe(summary), sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
