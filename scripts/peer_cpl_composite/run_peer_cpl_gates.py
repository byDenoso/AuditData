#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable

import numpy as np

import camb
from camb import dark_energy


@dataclass
class RunResult:
    label: str
    derived: dict[str, float]
    cls_path: str
    finite: bool


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def make_pars(model, lmax: int = 900):
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=70.0,
        ombh2=0.02237,
        omch2=0.1200,
        mnu=0.06,
        omk=0.0,
        tau=0.054,
        num_massive_neutrinos=1,
        nnu=3.046,
    )
    pars.InitPower.set_params(As=2.1e-9, ns=0.97)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    pars.WantTransfer = True
    pars.DarkEnergy = model
    return pars


def run_one(label: str, model, out: Path, lmax: int = 900) -> tuple[RunResult, np.ndarray]:
    pars = make_pars(model, lmax=lmax)
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)
    total = np.asarray(powers["total"], dtype=float)
    derived_raw = results.get_derived_params()
    derived = {
        key: float(derived_raw[key])
        for key in ("rdrag", "rstar", "zdrag", "zstar", "thetastar")
        if key in derived_raw
    }
    cls_path = out / f"{label}_cls.npy"
    np.save(cls_path, total)
    finite = bool(np.all(np.isfinite(total)) and all(np.isfinite(v) for v in derived.values()))
    return RunResult(label=label, derived=derived, cls_path=cls_path.name, finite=finite), total


def relative_spectrum_error(a: np.ndarray, b: np.ndarray, ell_min: int = 2) -> dict[str, float]:
    n = min(len(a), len(b))
    aa = a[ell_min:n]
    bb = b[ell_min:n]
    names = ("TT", "EE", "BB", "TE")
    out = {}
    for ix, name in enumerate(names):
        scale = max(float(np.max(np.abs(aa[:, ix]))), 1e-30)
        out[name] = float(np.max(np.abs(aa[:, ix] - bb[:, ix])) / scale)
    out["all_max"] = max(out.values())
    return out


def derived_error(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    out = {}
    for key in sorted(set(a) & set(b)):
        denom = max(abs(a[key]), 1e-30)
        out[key] = abs(a[key] - b[key]) / denom
    out["all_max"] = max(out.values(), default=0.0)
    return out


def stock_lambda():
    model = dark_energy.DarkEnergyPPF()
    model.set_params(w=-1.0, wa=0.0)
    return model


def composite_lambda():
    return dark_energy.PeerCPL().set_params(
        fde_zc=0.0,
        peer_enabled=False,
        cpl_enabled=True,
        w=-1.0,
        wa=0.0,
    )


def stock_cpl():
    model = dark_energy.DarkEnergyPPF()
    model.set_params(w=-0.9, wa=-0.2)
    return model


def composite_cpl():
    return dark_energy.PeerCPL().set_params(
        fde_zc=0.0,
        peer_enabled=False,
        cpl_enabled=True,
        w=-0.9,
        wa=-0.2,
    )


def stock_peer():
    model = dark_energy.EarlyQuintessence()
    model.set_params(
        n=3.0,
        theta_i=2.89155,
        use_zc=True,
        zc=10**3.81,
        fde_zc=0.088,
    )
    model.frac_lambda0 = 1.0
    return model


def composite_peer():
    return dark_energy.PeerCPL().set_params(
        n=3.0,
        theta_i=2.89155,
        zc=10**3.81,
        fde_zc=0.088,
        peer_enabled=True,
        cpl_enabled=True,
        w=-1.0,
        wa=0.0,
    )


def composite_active():
    return dark_energy.PeerCPL().set_params(
        n=3.0,
        theta_i=2.89155,
        zc=10**3.81,
        fde_zc=0.088,
        peer_enabled=True,
        cpl_enabled=True,
        w=-0.891075,
        wa=-0.227731,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lmax", type=int, default=900)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)

    factories: dict[str, Callable[[], object]] = {
        "stock_lambda": stock_lambda,
        "composite_lambda": composite_lambda,
        "stock_cpl": stock_cpl,
        "composite_cpl": composite_cpl,
        "stock_peer": stock_peer,
        "composite_peer": composite_peer,
        "composite_active": composite_active,
    }

    runs: dict[str, RunResult] = {}
    spectra: dict[str, np.ndarray] = {}
    for label, factory in factories.items():
        result, cls = run_one(label, factory(), out, lmax=args.lmax)
        runs[label] = result
        spectra[label] = cls

    comparisons = {
        "lambda_null": {
            "spectra": relative_spectrum_error(spectra["stock_lambda"], spectra["composite_lambda"]),
            "derived": derived_error(runs["stock_lambda"].derived, runs["composite_lambda"].derived),
            "threshold": 1e-8,
        },
        "cpl_null": {
            "spectra": relative_spectrum_error(spectra["stock_cpl"], spectra["composite_cpl"]),
            "derived": derived_error(runs["stock_cpl"].derived, runs["composite_cpl"].derived),
            "threshold": 2e-7,
        },
        "peer_null": {
            "spectra": relative_spectrum_error(spectra["stock_peer"], spectra["composite_peer"]),
            "derived": derived_error(runs["stock_peer"].derived, runs["composite_peer"].derived),
            "threshold": 5e-5,
        },
    }

    gates = {}
    for name, comp in comparisons.items():
        observed = max(comp["spectra"]["all_max"], comp["derived"]["all_max"])
        gates[name] = {
            "pass": bool(observed <= comp["threshold"]),
            "observed": observed,
            "threshold": comp["threshold"],
        }
    gates["active_finite"] = {
        "pass": runs["composite_active"].finite,
        "observed": runs["composite_active"].finite,
        "threshold": True,
    }

    payload = {
        "schema": "peer-cpl-composite-gate-1",
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "camb": camb.__version__,
            "lmax": args.lmax,
        },
        "runs": {key: asdict(value) for key, value in runs.items()},
        "comparisons": comparisons,
        "gates": gates,
        "all_pass": all(item["pass"] for item in gates.values()),
        "claim_boundary": (
            "Passing these gates validates composite CAMB wiring and null recovery. "
            "It is not a likelihood fit or model-selection result."
        ),
    }

    json_path = out / "peer_cpl_composite_gates.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# PEER + CPL composite CAMB gate",
        "",
        f"CAMB: `{camb.__version__}`",
        f"All gates pass: **{payload['all_pass']}**",
        "",
        "| Gate | Observed | Threshold | Pass |",
        "|---|---:|---:|---|",
    ]
    for name, gate in gates.items():
        report.append(f"| {name} | {gate['observed']} | {gate['threshold']} | {gate['pass']} |")
    report.extend([
        "",
        "## Active composite readout",
        "",
        "```json",
        json.dumps(runs["composite_active"].derived, indent=2, sort_keys=True),
        "```",
        "",
        "This report validates numerical wiring and null limits only.",
    ])
    report_path = out / "PEER_CPL_COMPOSITE_GATE.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    hashes = []
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{sha256(path)}  {path.name}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")

    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0 if payload["all_pass"] else 2)


if __name__ == "__main__":
    main()
