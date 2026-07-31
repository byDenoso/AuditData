from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np

BASE = {
    "H0": 70.795,
    "ombh2": 0.02290,
    "omch2": 0.12530,
    "tau": 0.0540,
    "As": 2.1086e-9,
    "ns": 0.9900,
    "mnu": 0.06,
    "nnu": 3.046,
    "Alens": 1.0,
    "log10_zc": 3.81,
    "theta_i": 2.89155,
}
STEPS = {
    "lnAs": 0.01,
    "ns": 0.005,
    "ombh2": 0.00020,
    "omch2": 0.0010,
    "tau": 0.005,
    "H0": 0.50,
    "Alens": 0.05,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def varied(anchor: dict, name: str, sign: int) -> dict:
    config = deepcopy(anchor)
    if name == "lnAs":
        config["As"] *= float(np.exp(sign * STEPS[name]))
    else:
        config[name] += sign * STEPS[name]
    return config


def make_params(config: dict, f_peer: float, n: float, lmax: int):
    import camb
    from camb.dark_energy import EarlyQuintessence

    if camb.__version__ != "1.6.6":
        raise RuntimeError(f"CAMB version mismatch: {camb.__version__}")
    pars = camb.CAMBparams()
    if f_peer > 0:
        dark_energy = EarlyQuintessence()
        dark_energy.set_params(
            n=n,
            theta_i=config["theta_i"],
            use_zc=True,
            zc=10.0 ** config["log10_zc"],
            fde_zc=f_peer,
        )
        pars.DarkEnergy = dark_energy
    pars.set_cosmology(
        H0=config["H0"],
        ombh2=config["ombh2"],
        omch2=config["omch2"],
        tau=config["tau"],
        mnu=config["mnu"],
        nnu=config["nnu"],
        num_massive_neutrinos=1,
    )
    pars.InitPower.set_params(As=config["As"], ns=config["ns"])
    pars.Alens = config["Alens"]
    pars.WantCls = True
    pars.Want_CMB = True
    pars.WantTransfer = False
    pars.NonLinear = camb.model.NonLinear_none
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    pars.DoLensing = True
    pars.Accuracy.AccuracyBoost = 1.2
    pars.Accuracy.lAccuracyBoost = 1.2
    pars.Accuracy.lSampleBoost = 1.2
    return pars


def run(output: Path, anchor_name: str, f_peer: float, n: float, lmax: int) -> None:
    import camb

    output.mkdir(parents=True, exist_ok=True)
    records = []
    for parameter in STEPS:
        for sign, suffix in ((-1, "minus"), (1, "plus")):
            config = varied(BASE, parameter, sign)
            results = camb.get_results(make_params(config, f_peer, n, lmax))
            lensed = np.asarray(
                results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True),
                dtype=float,
            )
            if lensed.shape[0] != lmax + 1 or lensed.shape[1] < 4:
                raise RuntimeError(f"unexpected spectrum shape {lensed.shape}")
            if not np.all(np.isfinite(lensed)):
                raise RuntimeError(f"non-finite derivative {parameter} {suffix}")
            label = f"deriv_{parameter}_{suffix}"
            npz = output / f"{label}.npz"
            meta = output / f"{label}.json"
            np.savez_compressed(npz, ell=np.arange(lmax + 1), lensed=lensed)
            metadata = {
                "test_id": "T-FN-001-TANGENT",
                "anchor_name": anchor_name,
                "anchor_f_peer": f_peer,
                "anchor_n": n,
                "parameter": parameter,
                "sign": sign,
                "step": STEPS[parameter],
                "config": config,
                "camb_version": camb.__version__,
                "lmax": lmax,
                "nonlinear_mode": "NonLinear_none",
            }
            meta.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
            records.append({
                "label": label,
                "npz_sha256": sha256(npz),
                "json_sha256": sha256(meta),
            })
            print("TANGENT_COMPLETE", anchor_name, label, flush=True)
    contract = {
        "test_id": "T-FN-001-TANGENT",
        "anchor_name": anchor_name,
        "anchor": {"f_peer": f_peer, "n": n},
        "canonical_config": {**BASE, "fde_zc": f_peer, "n": n},
        "derivative_steps": STEPS,
        "camb_version": camb.__version__,
        "derivative_spectra": len(records),
        "state": "TANGENT_BASIS_COMPLETED_AND_VALIDATED",
    }
    (output / "generation_contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--anchor-name", required=True)
    parser.add_argument("--f-peer", type=float, required=True)
    parser.add_argument("--n", type=float, required=True)
    parser.add_argument("--lmax", type=int, default=3000)
    args = parser.parse_args()
    run(args.output, args.anchor_name, args.f_peer, args.n, args.lmax)


if __name__ == "__main__":
    main()
