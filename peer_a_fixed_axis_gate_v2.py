#!/usr/bin/env python3
"""Planck PR3 fixed-axis gate for an exploratory PEER-A branch.

The full-mission SMICA-noSZ inpainted intensity field is used for the headline
map statistic. Half-mission intensity maps are used only as a stability check.
The gate calibrates a generic low-l anisotropy statistic with Gaussian
isotropic mocks. It does not implement an anisotropic PEER transfer function.
"""
from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import camb
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

FILES = {
    "full": "COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits",
    "hm1": "COM_CMB_IQU-smica-nosz_2048_R3.00_hm1.fits",
    "hm2": "COM_CMB_IQU-smica-nosz_2048_R3.00_hm2.fits",
}
SOURCES = {
    label: [
        f"https://irsawebops1.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/{name}",
        f"https://irsawebops2.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/{name}",
        f"https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/{name}",
    ]
    for label, name in FILES.items()
}
ACT_SPT_CONTEXT = {
    "source": "RESULTS_STATUS_20260731",
    "baseline": {"delta_chi2": -4.1868, "delta_aic": -2.1868},
    "plus_act": {"delta_chi2": -0.8576, "delta_aic": 1.1424},
    "plus_spt": {"delta_chi2": -2.7729, "delta_aic": -0.7729},
    "plus_act_spt": {"delta_chi2": 0.1574, "delta_aic": 2.1574},
    "boundary": "Additive stress profile; Planck/ACT/SPT cross-covariance not represented.",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)


def verify_inputs(data_dir: Path) -> list[dict[str, object]]:
    provenance: list[dict[str, object]] = []
    for label, filename in FILES.items():
        path = data_dir / filename
        if not path.is_file() or path.stat().st_size < 50_000_000:
            raise RuntimeError(f"Missing or truncated official input: {path}")
        provenance.append({
            "label": label,
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "candidate_sources": SOURCES[label],
        })
    return provenance


def angular_momentum_matrices(ell: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m_values = np.arange(-ell, ell + 1, dtype=float)
    raising = np.zeros((2 * ell + 1, 2 * ell + 1), dtype=complex)
    for index, m_value in enumerate(m_values[:-1]):
        raising[index + 1, index] = math.sqrt(ell * (ell + 1) - m_value * (m_value + 1))
    lowering = raising.conj().T
    return 0.5 * (raising + lowering), (raising - lowering) / (2j), np.diag(m_values)


def full_alm_vector(alm: np.ndarray, ell: int, lmax: int) -> np.ndarray:
    values = np.empty(2 * ell + 1, dtype=complex)
    for offset, m_value in enumerate(range(-ell, ell + 1)):
        if m_value >= 0:
            values[offset] = alm[hp.Alm.getidx(lmax, ell, m_value)]
        else:
            positive_m = -m_value
            positive = alm[hp.Alm.getidx(lmax, ell, positive_m)]
            values[offset] = ((-1) ** positive_m) * np.conj(positive)
    return values


def power_tensor_axis(alm: np.ndarray, ell: int, lmax: int) -> dict[str, object]:
    vector = full_alm_vector(alm, ell, lmax)
    norm = float(np.vdot(vector, vector).real)
    if not np.isfinite(norm) or norm <= 0:
        raise RuntimeError(f"Invalid multipole power at ell={ell}")
    operators = angular_momentum_matrices(ell)
    tensor = np.empty((3, 3), dtype=float)
    for i, first in enumerate(operators):
        for j, second in enumerate(operators):
            symmetric = 0.5 * (first @ second + second @ first)
            tensor[i, j] = float(np.vdot(vector, symmetric @ vector).real / norm)
    tensor = 0.5 * (tensor + tensor.T)
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    best = int(np.argmax(eigenvalues))
    axis = np.asarray(eigenvectors[:, best], dtype=float)
    axis /= np.linalg.norm(axis)
    return {
        "axis": axis,
        "concentration": float(eigenvalues[best] / np.trace(tensor)),
        "eigenvalues": eigenvalues,
    }


def undirected_angle(first: np.ndarray, second: np.ndarray) -> float:
    a = np.asarray(first, dtype=float) / np.linalg.norm(first)
    b = np.asarray(second, dtype=float) / np.linalg.norm(second)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a, b)), 0.0, 1.0))))


def combined_axis(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    a = np.asarray(first, dtype=float) / np.linalg.norm(first)
    b = np.asarray(second, dtype=float) / np.linalg.norm(second)
    if np.dot(a, b) < 0:
        b = -b
    output = a + b
    output /= np.linalg.norm(output)
    return output


def galactic(axis: np.ndarray) -> dict[str, float]:
    unit = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    return {
        "l_deg": float(np.degrees(np.arctan2(unit[1], unit[0])) % 360.0),
        "b_deg": float(np.degrees(np.arcsin(np.clip(unit[2], -1.0, 1.0)))),
    }


def fill_invalid(values: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values) & (np.abs(values) < 1.0e20) & (values != hp.UNSEEN)
    fraction_invalid = 1.0 - float(np.mean(valid))
    if not np.any(valid):
        raise RuntimeError("Map has no finite HEALPix samples")
    fill = float(np.median(values[valid]))
    output = values.copy()
    output[~valid] = fill
    return output, fraction_invalid


def load_alm(path: Path, field: int, lmax: int, nside_analysis: int) -> tuple[np.ndarray, float]:
    raw = hp.read_map(path, field=field, dtype=np.float64, memmap=True)
    cleaned, fraction_invalid = fill_invalid(raw)
    del raw
    cleaned = hp.ud_grade(cleaned, nside_out=nside_analysis, power=0)
    cleaned = hp.remove_dipole(cleaned, fitval=False)
    alm = hp.map2alm(cleaned, lmax=lmax, iter=3, use_pixel_weights=False)
    del cleaned
    gc.collect()
    return alm, fraction_invalid


def band_asymmetry(
    alm: np.ndarray,
    axis: np.ndarray,
    ell_min: int,
    ell_max: int,
    lmax: int,
    nside: int,
    pixel_vectors: np.ndarray,
) -> dict[str, float]:
    window = np.zeros(lmax + 1)
    window[ell_min : ell_max + 1] = 1.0
    band_map = hp.alm2map(hp.almxfl(alm, window), nside=nside, lmax=lmax)
    projection = axis @ pixel_vectors
    north = projection >= 0
    south = ~north
    north_power = float(np.mean(band_map[north] ** 2))
    south_power = float(np.mean(band_map[south] ** 2))
    amplitude = (north_power - south_power) / (north_power + south_power)
    return {
        "ell_min": ell_min,
        "ell_max": ell_max,
        "amplitude": float(amplitude),
        "abs_amplitude": float(abs(amplitude)),
        "north_power": north_power,
        "south_power": south_power,
    }


def analyze(alm: np.ndarray, lmax: int, nside_band: int, pixel_vectors: np.ndarray) -> dict[str, object]:
    quadrupole = power_tensor_axis(alm, 2, lmax)
    octopole = power_tensor_axis(alm, 3, lmax)
    axis = combined_axis(quadrupole["axis"], octopole["axis"])
    return {
        "quadrupole": {
            "axis": [float(value) for value in quadrupole["axis"]],
            "galactic": galactic(quadrupole["axis"]),
            "concentration": quadrupole["concentration"],
        },
        "octopole": {
            "axis": [float(value) for value in octopole["axis"]],
            "galactic": galactic(octopole["axis"]),
            "concentration": octopole["concentration"],
        },
        "alignment_angle_deg": undirected_angle(quadrupole["axis"], octopole["axis"]),
        "combined_axis": [float(value) for value in axis],
        "combined_galactic": galactic(axis),
        "low_band": band_asymmetry(alm, axis, 4, 30, lmax, nside_band, pixel_vectors),
        "control_band": band_asymmetry(alm, axis, 31, lmax, lmax, nside_band, pixel_vectors),
    }


def theory_cl(lmax: int) -> np.ndarray:
    parameters = camb.CAMBparams()
    parameters.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.120, tau=0.054)
    parameters.InitPower.set_params(As=2.1e-9, ns=0.965)
    parameters.set_for_lmax(max(lmax + 50, 120), lens_potential_accuracy=0)
    spectra = camb.get_results(parameters).get_cmb_power_spectra(parameters, CMB_unit="muK", raw_cl=True)
    cl = np.asarray(spectra["total"][: lmax + 1, 0], dtype=float)
    cl[:2] = 0.0
    if np.any(~np.isfinite(cl)) or np.any(cl[2:] <= 0):
        raise RuntimeError("Invalid CAMB TT spectrum")
    return cl


def empirical_p(values: np.ndarray, observed: float, lower: bool) -> float:
    count = np.count_nonzero(values <= observed) if lower else np.count_nonzero(values >= observed)
    return float((count + 1) / (values.size + 1))


def run_mocks(
    cl: np.ndarray,
    observed: dict[str, object],
    n_mock: int,
    seed: int,
    lmax: int,
    nside_band: int,
    pixel_vectors: np.ndarray,
    output_csv: Path,
) -> tuple[dict[str, object], list[dict[str, float]]]:
    np.random.seed(seed)
    rows: list[dict[str, float]] = []
    for index in range(n_mock):
        result = analyze(hp.synalm(cl, lmax=lmax, new=True), lmax, nside_band, pixel_vectors)
        rows.append({
            "mock": index,
            "alignment_angle_deg": float(result["alignment_angle_deg"]),
            "low_abs_amplitude": float(result["low_band"]["abs_amplitude"]),
            "control_abs_amplitude": float(result["control_band"]["abs_amplitude"]),
        })
        if (index + 1) % 128 == 0 or index + 1 == n_mock:
            with output_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            print(f"MOCK_PROGRESS {index + 1}/{n_mock}", flush=True)

    angles = np.array([row["alignment_angle_deg"] for row in rows])
    low = np.array([row["low_abs_amplitude"] for row in rows])
    control = np.array([row["control_abs_amplitude"] for row in rows])
    observed_angle = float(observed["alignment_angle_deg"])
    observed_low = float(observed["low_band"]["abs_amplitude"])
    observed_control = float(observed["control_band"]["abs_amplitude"])
    joint_count = np.count_nonzero((angles <= observed_angle) & (low >= observed_low))
    summary = {
        "n_mock": n_mock,
        "seed": seed,
        "p_alignment": empirical_p(angles, observed_angle, True),
        "p_low_modulation": empirical_p(low, observed_low, False),
        "p_control_modulation": empirical_p(control, observed_control, False),
        "p_joint": float((joint_count + 1) / (n_mock + 1)),
        "alignment_quantiles_deg": [float(value) for value in np.quantile(angles, [0.025, 0.16, 0.5, 0.84, 0.975])],
        "low_quantiles": [float(value) for value in np.quantile(low, [0.025, 0.16, 0.5, 0.84, 0.975])],
        "control_quantiles": [float(value) for value in np.quantile(control, [0.025, 0.16, 0.5, 0.84, 0.975])],
    }
    return summary, rows


def consistency(maps: dict[str, dict[str, object]]) -> dict[str, object]:
    full_axis = np.asarray(maps["full"]["combined_axis"])
    return {
        label: {
            "axis_angle_to_full_deg": undirected_angle(full_axis, np.asarray(maps[label]["combined_axis"])),
            "low_abs_difference": abs(float(maps[label]["low_band"]["abs_amplitude"]) - float(maps["full"]["low_band"]["abs_amplitude"])),
            "control_abs_difference": abs(float(maps[label]["control_band"]["abs_amplitude"]) - float(maps["full"]["control_band"]["abs_amplitude"])),
        }
        for label in ("hm1", "hm2")
    }


def decision(payload: dict[str, object]) -> dict[str, object]:
    full = payload["maps"]["full"]
    mocks = payload["mocks"]
    splits = payload["split_consistency"]
    stable_axis = all(float(splits[label]["axis_angle_to_full_deg"]) <= 15.0 for label in ("hm1", "hm2"))
    stable_low = all(float(splits[label]["low_abs_difference"]) <= 0.04 for label in ("hm1", "hm2"))
    anomalous = float(mocks["p_joint"]) < 0.05
    low = float(full["low_band"]["abs_amplitude"])
    control = float(full["control_band"]["abs_amplitude"])
    localized = control < max(0.02, 0.60 * low)
    if anomalous and stable_axis and stable_low and localized:
        status = "GENERIC_ANOMALY_ONLY"
    elif anomalous and stable_axis:
        status = "ANOMALY_NOT_LOW_L_LOCALIZED"
    else:
        status = "NO_PEER_A_SUPPORT_FROM_THIS_GATE"
    return {
        "status": status,
        "flags": {
            "joint_anomaly": anomalous,
            "half_mission_axis_stable": stable_axis,
            "half_mission_low_band_stable": stable_low,
            "low_l_localized": localized,
        },
        "claim_boundary": "A positive status remains a generic anisotropy result, not evidence that PEER caused the Axis of Evil.",
    }


def plot_results(output: Path, payload: dict[str, object], rows: list[dict[str, float]]) -> None:
    full = payload["maps"]["full"]
    angles = np.array([row["alignment_angle_deg"] for row in rows])
    low = np.array([row["low_abs_amplitude"] for row in rows])

    figure, axis = plt.subplots(figsize=(7, 4.2))
    axis.hist(angles, bins=50, density=True)
    axis.axvline(full["alignment_angle_deg"], linewidth=2, label=f"Planck {full['alignment_angle_deg']:.2f} deg")
    axis.set(xlabel="Quadrupole-octopole angle [deg]", ylabel="Mock density", title="Isotropic alignment calibration")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "alignment_calibration.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7, 4.2))
    axis.hist(low, bins=50, density=True)
    axis.axvline(full["low_band"]["abs_amplitude"], linewidth=2, label=f"Planck {full['low_band']['abs_amplitude']:.3f}")
    axis.set(xlabel="|hemispherical asymmetry|, ell=4-30", ylabel="Mock density", title="Fixed-axis low-l calibration")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "low_band_calibration.png", dpi=180)
    plt.close(figure)

    labels = ["full", "hm1", "hm2"]
    positions = np.arange(3)
    width = 0.36
    figure, axis = plt.subplots(figsize=(7, 4.2))
    axis.bar(positions - width / 2, [payload["maps"][label]["low_band"]["abs_amplitude"] for label in labels], width, label="ell=4-30")
    axis.bar(positions + width / 2, [payload["maps"][label]["control_band"]["abs_amplitude"] for label in labels], width, label="ell=31-64")
    axis.set_xticks(positions, labels)
    axis.set(ylabel="Absolute hemispherical asymmetry", title="Scale and half-mission consistency")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "split_consistency.png", dpi=180)
    plt.close(figure)


def write_report(output: Path, payload: dict[str, object]) -> None:
    full = payload["maps"]["full"]
    mocks = payload["mocks"]
    status = payload["decision"]["status"]
    report = f"""# PEER-A fixed-axis Planck gate

## 1. Estado observado

Planck PR3 SMICA-noSZ inpainted full-mission intensity gives a quadrupole-octopole power-tensor angle of **{full['alignment_angle_deg']:.4f} deg**. The combined undirected axis is **l={full['combined_galactic']['l_deg']:.3f} deg, b={full['combined_galactic']['b_deg']:.3f} deg**.

The fixed-axis asymmetry is **|A|={full['low_band']['abs_amplitude']:.6f}** for ell=4-30 and **|A|={full['control_band']['abs_amplitude']:.6f}** for ell=31-64.

## 2. Evidencia

- Gaussian isotropic mocks: **{mocks['n_mock']}**, seed **{mocks['seed']}**.
- Alignment empirical p: **{mocks['p_alignment']:.6g}**.
- Low-band modulation empirical p: **{mocks['p_low_modulation']:.6g}**.
- Joint empirical p: **{mocks['p_joint']:.6g}**.
- Full mission and both half missions use the same estimator. Half-mission invalid pixels are median-filled only for the stability diagnostic.

## 3. Inconsistencias

This gate uses one component-separation product, Gaussian isotropic mocks, temperature only, and an inpainted full-sky statistic. It does not include Planck FFP end-to-end simulations, alternative cleaners, polarization, or an anisotropic PEER transfer function. ACT/SPT enter only through the completed additive high-l profile, without cross-covariance.

## 4. Decisao

**{status}**

{payload['decision']['claim_boundary']}

## 5. Teste prioritario

Compute directional PEER transfer derivatives, project them into BipoSH coefficients, and compare PEER-A against a generic dipole-modulation model in Planck TT/TE/EE. ACT DR6 and SPT-3G then become high-l nulls with real window functions.

## 6. Impacto nos papers

No dipole or Axis-of-Evil claim enters the isotropic PEER paper. PEER-A remains a separate exploratory branch. A generic anomaly can motivate that branch, but cannot establish a new phenomenon.

## 7. Arquivos a atualizar

After artifact review: `ops/decision_log.md`, `ops/evidence_index.md`, `ops/run_registry.csv`, `papers/publication_matrix.md`, and only then `ops/peer_status.yaml` if promoted.
"""
    (output / "REPORT.md").write_text(report, encoding="utf-8")


def self_test() -> None:
    lmax = 8
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
    alm[hp.Alm.getidx(lmax, 3, 3)] = 1.0
    result = power_tensor_axis(alm, 3, lmax)
    assert abs(float(np.dot(result["axis"], [0.0, 0.0, 1.0]))) > 0.999999
    assert undirected_angle(np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])) < 1.0e-10
    print("SELF_TEST_PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/peer_a_fixed_axis"))
    parser.add_argument("--data-dir", type=Path, default=Path("data/planck_pr3_smica_nosz"))
    parser.add_argument("--n-mock", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--lmax", type=int, default=64)
    parser.add_argument("--nside-analysis", type=int, default=128)
    parser.add_argument("--nside-band", type=int, default=32)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.lmax < 32 or args.nside_band * 3 - 1 < args.lmax:
        raise SystemExit("Invalid lmax/nside contract")

    started = time.time()
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = verify_inputs(args.data_dir)
    pixel_vectors = np.array(hp.pix2vec(args.nside_band, np.arange(hp.nside2npix(args.nside_band))))

    maps: dict[str, dict[str, object]] = {}
    invalid_fractions: dict[str, float] = {}
    fields = {"full": 1, "hm1": 0, "hm2": 0}
    for label in ("full", "hm1", "hm2"):
        print(f"MAP_ANALYSIS {label}", flush=True)
        alm, invalid_fraction = load_alm(args.data_dir / FILES[label], fields[label], args.lmax, args.nside_analysis)
        maps[label] = analyze(alm, args.lmax, args.nside_band, pixel_vectors)
        invalid_fractions[label] = invalid_fraction
        del alm
        gc.collect()

    mock_summary, mock_rows = run_mocks(
        theory_cl(args.lmax), maps["full"], args.n_mock, args.seed,
        args.lmax, args.nside_band, pixel_vectors, args.output / "mocks.csv",
    )
    payload: dict[str, object] = {
        "evidence_class": "PLANCK_PR3_SMICA_NOSZ_MAP_GAUSSIAN_NULL_GATE",
        "claim_boundary": "Generic statistical-anisotropy gate; not an anisotropic PEER likelihood or a new-phenomenon detection.",
        "configuration": {
            "n_mock": args.n_mock,
            "seed": args.seed,
            "lmax": args.lmax,
            "nside_analysis": args.nside_analysis,
            "nside_band": args.nside_band,
            "full_field": "I_STOKES_INP (field 1)",
            "half_mission_field": "I_STOKES (field 0)",
        },
        "maps": maps,
        "invalid_pixel_fraction_before_fill": invalid_fractions,
        "mocks": mock_summary,
        "split_consistency": consistency(maps),
        "act_spt_context": ACT_SPT_CONTEXT,
        "provenance": provenance,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "healpy": hp.__version__,
            "camb": camb.__version__,
            "github_sha": os.environ.get("GITHUB_SHA"),
            "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        },
    }
    payload["decision"] = decision(payload)
    payload["elapsed_seconds"] = time.time() - started
    write_json(args.output / "results.json", payload)
    write_json(args.output / "provenance.json", provenance)
    plot_results(args.output, payload, mock_rows)
    write_report(args.output, payload)
    summary = {
        "decision": payload["decision"]["status"],
        "axis": maps["full"]["combined_galactic"],
        "angle_deg": maps["full"]["alignment_angle_deg"],
        "low_abs": maps["full"]["low_band"]["abs_amplitude"],
        "control_abs": maps["full"]["control_band"]["abs_amplitude"],
        "p_alignment": mock_summary["p_alignment"],
        "p_low": mock_summary["p_low_modulation"],
        "p_joint": mock_summary["p_joint"],
        "splits": payload["split_consistency"],
    }
    print("SCIENTIFIC_SUMMARY=" + json.dumps(summary, sort_keys=True), flush=True)
    print("GATE_COMPLETE", payload["decision"]["status"], flush=True)


if __name__ == "__main__":
    main()
