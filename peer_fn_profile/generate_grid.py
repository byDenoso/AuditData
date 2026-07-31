from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

CANONICAL = {
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


def label_for(f_peer: float, n: float) -> str:
    return f"f{int(round(f_peer * 1000)):03d}_n{int(round(n * 100)):03d}"


def parse_grid(text: str) -> list[float]:
    values = [float(token.strip()) for token in text.split(",") if token.strip()]
    array = np.asarray(values, dtype=float)
    if not values or not np.all(np.isfinite(array)):
        raise ValueError("f grid must be non-empty and finite")
    if np.any(np.diff(array) <= 0):
        raise ValueError("f grid must be strictly increasing")
    if values[0] < 0 or values[-1] > 0.14:
        raise ValueError("f grid must stay inside [0, 0.14]")
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def make_params(f_peer: float, n: float, lmax: int):
    import camb
    from camb.dark_energy import EarlyQuintessence

    if camb.__version__ != "1.6.6":
        raise RuntimeError(f"CAMB version mismatch: {camb.__version__}")
    pars = camb.CAMBparams()
    if f_peer > 0:
        dark_energy = EarlyQuintessence()
        dark_energy.set_params(
            n=n,
            theta_i=CANONICAL["theta_i"],
            use_zc=True,
            zc=10.0 ** CANONICAL["log10_zc"],
            fde_zc=f_peer,
        )
        pars.DarkEnergy = dark_energy
    pars.set_cosmology(
        H0=CANONICAL["H0"],
        ombh2=CANONICAL["ombh2"],
        omch2=CANONICAL["omch2"],
        tau=CANONICAL["tau"],
        mnu=CANONICAL["mnu"],
        nnu=CANONICAL["nnu"],
        num_massive_neutrinos=1,
    )
    pars.InitPower.set_params(As=CANONICAL["As"], ns=CANONICAL["ns"])
    pars.Alens = CANONICAL["Alens"]
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


def run_point(output: Path, f_peer: float, n: float, lmax: int) -> dict:
    import camb

    label = label_for(f_peer, n)
    results = camb.get_results(make_params(f_peer, n, lmax))
    lensed = np.asarray(
        results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True),
        dtype=float,
    )
    ell = np.arange(lmax + 1, dtype=int)
    if lensed.shape[0] != lmax + 1 or lensed.shape[1] < 4:
        raise RuntimeError(f"unexpected CAMB spectrum shape {lensed.shape}")
    if not np.all(np.isfinite(lensed)):
        raise RuntimeError(f"non-finite spectrum for {label}")
    derived = results.get_derived_params()
    metadata = {
        "test_id": "T-FN-001",
        "label": label,
        "model": "lcdm_null_boundary" if f_peer == 0 else "early_quintessence",
        "n_identified": f_peer > 0,
        "camb_version": camb.__version__,
        "config": {**CANONICAL, "fde_zc": f_peer, "n": n},
        "zc": 10.0 ** CANONICAL["log10_zc"],
        "lmax": lmax,
        "nonlinear_mode": "NonLinear_none",
        "do_lensing": True,
        "raw_cl": True,
        "cmb_unit": "muK",
        "derived": {
            key: float(derived[key])
            for key in ("zstar", "rstar", "thetastar", "zdrag", "rdrag", "zeq")
            if key in derived
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    npz_path = output / f"{label}.npz"
    json_path = output / f"{label}.json"
    np.savez_compressed(npz_path, ell=ell, lensed=lensed)
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    print(
        "GRID_POINT_COMPLETE",
        json.dumps({"label": label, "rdrag": metadata["derived"].get("rdrag")}, sort_keys=True),
        flush=True,
    )
    return {**metadata, "npz_sha256": sha256(npz_path), "json_sha256": sha256(json_path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=float, required=True)
    parser.add_argument("--f-grid", required=True)
    parser.add_argument("--lmax", type=int, default=3000)
    args = parser.parse_args()
    if not 2.8 <= args.n <= 4.5:
        raise SystemExit("n must be inside [2.8, 4.5]")
    if args.lmax < 2508:
        raise SystemExit("lmax must be >= 2508")
    f_values = parse_grid(args.f_grid)
    points = [run_point(args.output, f_peer, args.n, args.lmax) for f_peer in f_values]
    code = int(round(args.n * 100))
    contract = {
        "test_id": "T-FN-001",
        "row_n": args.n,
        "f_values": f_values,
        "point_count": len(points),
        "canonical": CANONICAL,
        "camb_version": points[0]["camb_version"],
        "completion": "ROW_COMPLETE",
    }
    (args.output / f"row_n{code:03d}_contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
