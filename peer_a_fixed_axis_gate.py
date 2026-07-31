#!/usr/bin/env python3
"""PEER-A fixed-axis Planck map gate.

Evidence class: Planck PR3 map-level diagnostic with Gaussian isotropic mocks.
Claim boundary: this tests generic low-l statistical anisotropy. It does not
implement an anisotropic PEER Boltzmann hierarchy and cannot identify PEER as
the cause of an anomaly.
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
import urllib.request

import camb
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

PLANCK_BASE = "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb"
PLANCK_FILES = {
    "full": "COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits",
    "hm1": "COM_CMB_IQU-smica-nosz_2048_R3.00_hm1.fits",
    "hm2": "COM_CMB_IQU-smica-nosz_2048_R3.00_hm2.fits",
}
ACT_SPT_CONTEXT = {
    "source": "RESULTS_STATUS_20260731",
    "baseline": {"delta_chi2": -4.1868, "delta_aic": -2.1868},
    "plus_act": {"delta_chi2": -0.8576, "delta_aic": 1.1424},
    "plus_spt": {"delta_chi2": -2.7729, "delta_aic": -0.7729},
    "plus_act_spt": {"delta_chi2": 0.1574, "delta_aic": 2.1574},
    "boundary": "Additive stress profile; Planck/ACT/SPT cross-covariance not represented.",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)


def download(url: str, target: Path, minimum_bytes: int = 50_000_000) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size >= minimum_bytes:
        return {
            "url": url,
            "path": str(target),
            "bytes": target.stat().st_size,
            "sha256": sha256_file(target),
            "cache_hit": True,
        }
    partial = target.with_suffix(target.suffix + ".part")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "PEER-A-gate/1.0"})
    with urllib.request.urlopen(request, timeout=240) as response, partial.open("wb") as output:
        while block := response.read(16 * 1024 * 1024):
            output.write(block)
    size = partial.stat().st_size
    if size < minimum_bytes:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Download truncated: {url}, {size} bytes")
    partial.replace(target)
    return {
        "url": url,
        "path": str(target),
        "bytes": size,
        "sha256": sha256_file(target),
        "cache_hit": False,
    }


def angular_momentum_matrices(ell: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m_values = np.arange(-ell, ell + 1, dtype=float)
    size = 2 * ell + 1
    raising = np.zeros((size, size), dtype=complex)
    for index, m_value in enumerate(m_values[:-1]):
        raising[index + 1, index] = math.sqrt(ell * (ell + 1) - m_value * (m_value + 1))
    lowering = raising.conj().T
    return 0.5 * (raising + lowering), (raising - lowering) / (2j), np.diag(m_values)


def full_alm_vector(alm: np.ndarray, ell: int, lmax: int) -> np.ndarray:
    vector = np.empty(2 * ell + 1, dtype=complex)
    for offset, m_value in enumerate(range(-ell, ell + 1)):
        if m_value >= 0:
            vector[offset] = alm[hp.Alm.getidx(lmax, ell, m_value)]
        else:
            positive_m = -m_value
            positive = alm[hp.Alm.getidx(lmax, ell, positive_m)]
            vector[offset] = ((-1) ** positive_m) * np.conj(positive)
    return vector


def power_tensor_axis(alm: np.ndarray, ell: int, lmax: int) -> dict[str, object]:
    vector = full_alm_vector(alm, ell, lmax)
    norm = float(np.vdot(vector, vector).real)
    if not np.isfinite(norm) or norm <= 0:
        raise RuntimeError(f"Invalid alm power at ell={ell}")
    matrices = angular_momentum_matrices(ell)
    tensor = np.empty((3, 3), dtype=float)
    for i, first in enumerate(matrices):
        for j, second in enumerate(matrices):
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


def undirected_angle_deg(first: np.ndarray, second: np.ndarray) -> float:
    a = np.asarray(first, dtype=float) / np.linalg.norm(first)
    b = np.asarray(second, dtype=float) / np.linalg.norm(second)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a, b)), 0.0, 1.0))))


def combine_axes(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    a = np.asarray(first, dtype=float) / np.linalg.norm(first)
    b = np.asarray(second, dtype=float) / np.linalg.norm(second)
    if np.dot(a, b) < 0:
        b = -b
    combined = a + b
    combined /= np.linalg.norm(combined)
    return combined


def galactic_coordinates(axis: np.ndarray) -> dict[str, float]:
    unit = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    return {
        "l_deg": float(np.degrees(np.arctan2(unit[1], unit[0])) % 360.0),
        "b_deg": float(np.degrees(np.arcsin(np.clip(unit[2], -1.0, 1.0)))),
    }


def band_asymmetry(
    alm: np.ndarray,
    axis: np.ndarray,
    ell_min: int,
    ell_max: int,
    lmax: int,
    nside: int,
    mask: np.ndarray,
    pixel_vectors: np.ndarray,
) -> dict[str, float]:
    window = np.zeros(lmax + 1)
    window[ell_min : ell_max + 1] = 1.0
    filtered_alm = hp.almxfl(alm, window, inplace=False)
    filtered_map = hp.alm2map(filtered_alm, nside=nside, lmax=lmax)
    projection = axis @ pixel_vectors
    valid = mask & np.isfinite(filtered_map)
    north = valid & (projection >= 0)
    south = valid & (projection < 0)
    if np.count_nonzero(north) < 100 or np.count_nonzero(south) < 100:
        raise RuntimeError("Insufficient unmasked pixels in one hemisphere")
    north_power = float(np.mean(filtered_map[north] ** 2))
    south_power = float(np.mean(filtered_map[south] ** 2))
    amplitude = (north_power - south_power) / (north_power + south_power)
    return {
        "ell_min": ell_min,
        "ell_max": ell_max,
        "amplitude": float(amplitude),
        "abs_amplitude": float(abs(amplitude)),
        "north_power": north_power,
        "south_power": south_power,
        "north_pixels": int(np.count_nonzero(north)),
        "south_pixels": int(np.count_nonzero(south)),
    }


def analyze_alm(
    alm: np.ndarray,
    lmax: int,
    nside_band: int,
    mask: np.ndarray,
    pixel_vectors: np.ndarray,
) -> dict[str, object]:
    quadrupole = power_tensor_axis(alm, 2, lmax)
    octopole = power_tensor_axis(alm, 3, lmax)
    combined = combine_axes(quadrupole["axis"], octopole["axis"])
    return {
        "quadrupole": {
            "axis": [float(value) for value in quadrupole["axis"]],
            "galactic": galactic_coordinates(quadrupole["axis"]),
            "concentration": quadrupole["concentration"],
        },
        "octopole": {
            "axis": [float(value) for value in octopole["axis"]],
            "galactic": galactic_coordinates(octopole["axis"]),
            "concentration": octopole["concentration"],
        },
        "alignment_angle_deg": undirected_angle_deg(quadrupole["axis"], octopole["axis"]),
        "combined_axis": [float(value) for value in combined],
        "combined_galactic": galactic_coordinates(combined),
        "low_band": band_asymmetry(alm, combined, 4, 30, lmax, nside_band, mask, pixel_vectors),
        "control_band": band_asymmetry(alm, combined, 31, lmax, lmax, nside_band, mask, pixel_vectors),
    }


def load_temperature_alm(path: Path, lmax: int, nside_analysis: int) -> np.ndarray:
    temperature = hp.read_map(path, field=0, dtype=np.float64, memmap=True)
    temperature = hp.ud_grade(temperature, nside_out=nside_analysis, power=0)
    temperature[~np.isfinite(temperature)] = hp.UNSEEN
    temperature = hp.remove_dipole(temperature, bad=hp.UNSEEN, fitval=False)
    alm = hp.map2alm(temperature, lmax=lmax, iter=3, use_pixel_weights=False)
    del temperature
    gc.collect()
    return alm


def load_mask(full_map_path: Path, nside: int) -> np.ndarray:
    mask = hp.read_map(full_map_path, field=3, dtype=np.float32, memmap=True)
    mask = hp.ud_grade(mask, nside_out=nside, power=0)
    return np.asarray(mask >= 0.90, dtype=bool)


def theory_tt_cl(lmax: int) -> np.ndarray:
    parameters = camb.CAMBparams()
    parameters.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.120, tau=0.054)
    parameters.InitPower.set_params(As=2.1e-9, ns=0.965)
    parameters.set_for_lmax(max(120, lmax + 50), lens_potential_accuracy=0)
    results = camb.get_results(parameters)
    spectra = results.get_cmb_power_spectra(parameters, CMB_unit="muK", raw_cl=True)
    cl = np.asarray(spectra["total"][: lmax + 1, 0], dtype=float)
    cl[:2] = 0.0
    if np.any(~np.isfinite(cl)) or np.any(cl[2:] <= 0):
        raise RuntimeError("CAMB produced an invalid TT spectrum")
    return cl


def empirical_p(values: np.ndarray, observed: float, lower_tail: bool) -> float:
    count = np.count_nonzero(values <= observed) if lower_tail else np.count_nonzero(values >= observed)
    return float((count + 1) / (values.size + 1))


def run_mocks(
    cl: np.ndarray,
    observed: dict[str, object],
    n_mock: int,
    seed: int,
    lmax: int,
    nside_band: int,
    mask: np.ndarray,
    pixel_vectors: np.ndarray,
    csv_path: Path,
) -> tuple[dict[str, object], list[dict[str, float]]]:
    np.random.seed(seed)
    rows: list[dict[str, float]] = []
    for index in range(n_mock):
        alm = hp.synalm(cl, lmax=lmax, new=True)
        result = analyze_alm(alm, lmax, nside_band, mask, pixel_vectors)
        rows.append({
            "mock": index,
            "alignment_angle_deg": float(result["alignment_angle_deg"]),
            "low_abs_amplitude": float(result["low_band"]["abs_amplitude"]),
            "control_abs_amplitude": float(result["control_band"]["abs_amplitude"]),
        })
        if (index + 1) % 128 == 0 or index + 1 == n_mock:
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
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
    joint = np.count_nonzero((angles <= observed_angle) & (low >= observed_low))
    summary = {
        "n_mock": n_mock,
        "seed": seed,
        "p_alignment": empirical_p(angles, observed_angle, lower_tail=True),
        "p_low_modulation": empirical_p(low, observed_low, lower_tail=False),
        "p_control_modulation": empirical_p(control, observed_control, lower_tail=False),
        "p_joint": float((joint + 1) / (n_mock + 1)),
        "alignment_quantiles_deg": [float(value) for value in np.quantile(angles, [0.025, 0.16, 0.5, 0.84, 0.975])],
        "low_quantiles": [float(value) for value in np.quantile(low, [0.025, 0.16, 0.5, 0.84, 0.975])],
        "control_quantiles": [float(value) for value in np.quantile(control, [0.025, 0.16, 0.5, 0.84, 0.975])],
    }
    return summary, rows


def split_consistency(maps: dict[str, dict[str, object]]) -> dict[str, object]:
    full_axis = np.asarray(maps["full"]["combined_axis"])
    output: dict[str, object] = {}
    for label in ("hm1", "hm2"):
        output[label] = {
            "axis_angle_to_full_deg": undirected_angle_deg(full_axis, np.asarray(maps[label]["combined_axis"])),
            "low_abs_amplitude_difference": abs(float(maps[label]["low_band"]["abs_amplitude"]) - float(maps["full"]["low_band"]["abs_amplitude"])),
            "control_abs_amplitude_difference": abs(float(maps[label]["control_band"]["abs_amplitude"]) - float(maps["full"]["control_band"]["abs_amplitude"])),
        }
    return output


def decide(payload: dict[str, object]) -> dict[str, object]:
    full = payload["maps"]["full"]
    mocks = payload["mocks"]
    splits = payload["split_consistency"]
    stable_axis = all(float(splits[label]["axis_angle_to_full_deg"]) <= 15.0 for label in ("hm1", "hm2"))
    stable_low = all(float(splits[label]["low_abs_amplitude_difference"]) <= 0.04 for label in ("hm1", "hm2"))
    anomaly = float(mocks["p_joint"]) < 0.05
    low = float(full["low_band"]["abs_amplitude"])
    control = float(full["control_band"]["abs_amplitude"])
    localized = control < max(0.02, 0.60 * low)
    if anomaly and stable_axis and stable_low and localized:
        status = "GENERIC_ANOMALY_ONLY"
    elif anomaly and stable_axis:
        status = "ANOMALY_NOT_LOW_L_LOCALIZED"
    else:
        status = "NO_PEER_A_SUPPORT_FROM_THIS_GATE"
    return {
        "status": status,
        "flags": {
            "joint_anomaly": anomaly,
            "half_mission_axis_stable": stable_axis,
            "half_mission_low_band_stable": stable_low,
            "low_l_localized": localized,
        },
        "claim_boundary": "Even a positive result is a generic map anomaly, not evidence that PEER caused the Axis of Evil.",
    }


def make_plots(output: Path, payload: dict[str, object], rows: list[dict[str, float]]) -> None:
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
    low_values = [payload["maps"][label]["low_band"]["abs_amplitude"] for label in labels]
    controls = [payload["maps"][label]["control_band"]["abs_amplitude"] for label in labels]
    positions = np.arange(3)
    width = 0.36
    figure, axis = plt.subplots(figsize=(7, 4.2))
    axis.bar(positions - width / 2, low_values, width, label="ell=4-30")
    axis.bar(positions + width / 2, controls, width, label="ell=31-64")
    axis.set_xticks(positions, labels)
    axis.set(ylabel="Absolute hemispherical asymmetry", title="Scale and half-mission consistency")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "split_consistency.png", dpi=180)
    plt.close(figure)


def write_report(output: Path, payload: dict[str, object]) -> None:
    full = payload["maps"]["full"]
    mocks = payload["mocks"]
    decision = payload["decision"]
    text = f"""# PEER-A fixed-axis Planck gate

## 1. Estado observado

Planck PR3 SMICA-noSZ gives a quadrupole-octopole power-tensor angle of **{full['alignment_angle_deg']:.4f} deg**. The combined undirected axis is **l={full['combined_galactic']['l_deg']:.3f} deg, b={full['combined_galactic']['b_deg']:.3f} deg**.

The fixed-axis power asymmetry is **|A|={full['low_band']['abs_amplitude']:.6f}** for ell=4-30 and **|A|={full['control_band']['abs_amplitude']:.6f}** for ell=31-64.

## 2. Evidencia

- Gaussian isotropic mocks: **{mocks['n_mock']}**.
- Alignment empirical p: **{mocks['p_alignment']:.6g}**.
- Low-band modulation empirical p: **{mocks['p_low_modulation']:.6g}**.
- Joint empirical p: **{mocks['p_joint']:.6g}**.
- Full mission and both half missions were processed with the same estimator and mask.

## 3. Inconsistencias

This is a Gaussian map-level gate. It does not include Planck FFP end-to-end simulations, alternative component-separation maps, polarization, or an anisotropic PEER transfer function. ACT/SPT enter only through the already completed additive high-l profile, whose cross-covariance is not represented.

## 4. Decisao

**{decision['status']}**

{decision['claim_boundary']}

## 5. Teste prioritario

Calculate directional PEER transfer derivatives, map them to BipoSH coefficients, and compare a PEER-specific response against a generic dipole-modulation model in Planck TT/TE/EE. ACT DR6 and SPT-3G then become map/likelihood high-l nulls.

## 6. Impacto nos papers

The current isotropic PEER paper receives no dipole or Axis-of-Evil claim. PEER-A remains a separate exploratory branch. A positive generic anomaly can justify the branch, but cannot identify a new phenomenon.

## 7. Arquivos a atualizar

After artifact review: `ops/decision_log.md`, `ops/evidence_index.md`, `ops/run_registry.csv`, `papers/publication_matrix.md`, and then `ops/peer_status.yaml` if canonically promoted.
"""
    (output / "REPORT.md").write_text(text, encoding="utf-8")


def self_test() -> None:
    lmax = 8
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
    alm[hp.Alm.getidx(lmax, 3, 3)] = 1.0
    result = power_tensor_axis(alm, 3, lmax)
    assert abs(float(np.dot(result["axis"], [0.0, 0.0, 1.0]))) > 0.999999
    assert undirected_angle_deg(np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])) < 1e-10
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
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test()
        return
    if arguments.lmax < 32 or arguments.nside_band * 3 - 1 < arguments.lmax:
        raise SystemExit("Invalid lmax/nside-band contract")

    started = time.time()
    arguments.output.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    provenance: list[dict[str, object]] = []
    for label, filename in PLANCK_FILES.items():
        target = arguments.data_dir / filename
        source = f"{PLANCK_BASE}/{filename}"
        provenance.append(download(source, target))
        paths[label] = target

    mask = load_mask(paths["full"], arguments.nside_band)
    pixel_vectors = np.array(hp.pix2vec(arguments.nside_band, np.arange(hp.nside2npix(arguments.nside_band))))
    maps: dict[str, dict[str, object]] = {}
    for label in ("full", "hm1", "hm2"):
        print(f"MAP_ANALYSIS {label}", flush=True)
        alm = load_temperature_alm(paths[label], arguments.lmax, arguments.nside_analysis)
        maps[label] = analyze_alm(alm, arguments.lmax, arguments.nside_band, mask, pixel_vectors)
        del alm
        gc.collect()

    cl = theory_tt_cl(arguments.lmax)
    mock_summary, mock_rows = run_mocks(
        cl=cl,
        observed=maps["full"],
        n_mock=arguments.n_mock,
        seed=arguments.seed,
        lmax=arguments.lmax,
        nside_band=arguments.nside_band,
        mask=mask,
        pixel_vectors=pixel_vectors,
        csv_path=arguments.output / "mocks.csv",
    )

    payload: dict[str, object] = {
        "evidence_class": "PLANCK_PR3_MAP_LEVEL_GAUSSIAN_NULL_GATE",
        "claim_boundary": "Generic statistical-anisotropy gate; not an anisotropic PEER likelihood or a new-phenomenon detection.",
        "configuration": {
            "n_mock": arguments.n_mock,
            "seed": arguments.seed,
            "lmax": arguments.lmax,
            "nside_analysis": arguments.nside_analysis,
            "nside_band": arguments.nside_band,
        },
        "maps": maps,
        "mocks": mock_summary,
        "split_consistency": split_consistency(maps),
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
    payload["decision"] = decide(payload)
    payload["elapsed_seconds"] = time.time() - started
    write_json(arguments.output / "results.json", payload)
    write_json(arguments.output / "provenance.json", provenance)
    make_plots(arguments.output, payload, mock_rows)
    write_report(arguments.output, payload)
    print("SCIENTIFIC_SUMMARY=" + json.dumps({
        "decision": payload["decision"]["status"],
        "axis": maps["full"]["combined_galactic"],
        "angle_deg": maps["full"]["alignment_angle_deg"],
        "low_abs": maps["full"]["low_band"]["abs_amplitude"],
        "control_abs": maps["full"]["control_band"]["abs_amplitude"],
        "p_alignment": mock_summary["p_alignment"],
        "p_low": mock_summary["p_low_modulation"],
        "p_joint": mock_summary["p_joint"],
    }, sort_keys=True), flush=True)
    print("GATE_COMPLETE", payload["decision"]["status"], flush=True)


if __name__ == "__main__":
    main()
