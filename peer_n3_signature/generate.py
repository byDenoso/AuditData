from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


@dataclass(frozen=True)
class CanonicalConfig:
    H0: float = 70.795
    ombh2: float = 0.02290
    omch2: float = 0.12530
    tau: float = 0.0540
    As: float = 2.1086e-9
    ns: float = 0.9900
    mnu: float = 0.06
    nnu: float = 3.046
    Alens: float = 1.0
    n: float = 3.0
    fde_zc: float = 0.0880
    log10_zc: float = 3.81
    theta_i: float = 2.89155

    @property
    def zc(self) -> float:
        return 10.0 ** self.log10_zc


DERIVATIVE_STEPS: dict[str, float] = {
    "lnAs": 0.01,
    "ns": 0.005,
    "ombh2": 0.00020,
    "omch2": 0.0010,
    "tau": 0.005,
    "H0": 0.50,
    "Alens": 0.05,
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


def import_camb():
    import camb
    from camb.dark_energy import AxionEffectiveFluid, EarlyQuintessence

    if camb.__version__ != "1.6.6":
        raise RuntimeError(f"CAMB version mismatch: {camb.__version__}")
    return camb, AxionEffectiveFluid, EarlyQuintessence


def make_params(
    config: CanonicalConfig,
    model: str,
    lmax: int,
    fluid_params: dict[str, float] | None = None,
    want_cls: bool = True,
):
    camb, AxionEffectiveFluid, EarlyQuintessence = import_camb()
    pars = camb.CAMBparams()
    if model == "scalar":
        dark_energy = EarlyQuintessence()
        dark_energy.set_params(
            n=config.n,
            theta_i=config.theta_i,
            use_zc=True,
            zc=config.zc,
            fde_zc=config.fde_zc,
        )
        pars.DarkEnergy = dark_energy
    elif model == "fluid":
        if fluid_params is None:
            raise ValueError("fluid_params are required for fluid model")
        dark_energy = AxionEffectiveFluid()
        dark_energy.set_params(
            w_n=0.5,
            theta_i=config.theta_i,
            zc=float(fluid_params["zc"]),
            fde_zc=float(fluid_params["fde_zc"]),
        )
        pars.DarkEnergy = dark_energy
    elif model != "lcdm":
        raise ValueError(f"unknown model: {model}")

    pars.set_cosmology(
        H0=config.H0,
        ombh2=config.ombh2,
        omch2=config.omch2,
        tau=config.tau,
        mnu=config.mnu,
        nnu=config.nnu,
        num_massive_neutrinos=1,
    )
    pars.InitPower.set_params(As=config.As, ns=config.ns)
    pars.Alens = config.Alens
    pars.WantCls = want_cls
    pars.Want_CMB = want_cls
    pars.WantTransfer = False
    pars.NonLinear = camb.model.NonLinear_none
    if want_cls:
        pars.set_for_lmax(lmax, lens_potential_accuracy=0)
        pars.Accuracy.AccuracyBoost = 1.2
        pars.Accuracy.lAccuracyBoost = 1.2
        pars.Accuracy.lSampleBoost = 1.2
    return pars


def background_from_results(results, z: np.ndarray) -> dict[str, Any]:
    derived = results.get_derived_params()
    selected = {
        key: float(derived[key])
        for key in ("zstar", "rstar", "thetastar", "zdrag", "rdrag", "kd", "zeq")
        if key in derived
    }
    selected["H0"] = float(results.Params.H0)
    return {
        "z": z,
        "H": np.asarray(results.hubble_parameter(z), dtype=float),
        "omega_de": np.asarray(results.get_Omega("de", z=z), dtype=float),
        "derived": selected,
    }


def run_spectrum(
    output_dir: Path,
    label: str,
    config: CanonicalConfig,
    model: str,
    lmax: int,
    fluid_params: dict[str, float] | None = None,
) -> dict[str, Any]:
    camb, _, _ = import_camb()
    pars = make_params(config, model=model, lmax=lmax, fluid_params=fluid_params, want_cls=True)
    results = camb.get_results(pars)
    lensed = np.asarray(results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True))
    unlensed = np.asarray(results.get_unlensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True))
    lens = np.asarray(results.get_lens_potential_cls(lmax=lmax, raw_cl=True))
    z_background = np.unique(np.concatenate(([0.0], np.logspace(-1, 5.2, 360))))
    background = background_from_results(results, z_background)
    arrays = {
        "ell": np.arange(lmax + 1, dtype=int),
        "lensed": lensed,
        "unlensed": unlensed,
        "lens": lens,
        "z": background["z"],
        "H": background["H"],
        "omega_de": background["omega_de"],
    }
    for key, array in arrays.items():
        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"non-finite {key} for {label}")
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / f"{label}.npz", **arrays)
    metadata = {
        "label": label,
        "model": model,
        "camb_version": camb.__version__,
        "config": asdict(config),
        "zc": config.zc,
        "fluid_params": fluid_params,
        "derived": background["derived"],
        "nonlinear_mode": "NonLinear_none",
    }
    (output_dir / f"{label}.json").write_text(json.dumps(json_safe(metadata), indent=2, sort_keys=True), encoding="utf-8")
    print("MODEL_COMPLETE", label, json.dumps(metadata["derived"], sort_keys=True), flush=True)
    return metadata


def fit_effective_fluid(config: CanonicalConfig, target_npz: Path, lmax: int) -> dict[str, float]:
    camb, _, _ = import_camb()
    target = np.load(target_npz)
    z_all = np.asarray(target["z"], dtype=float)
    omega_all = np.asarray(target["omega_de"], dtype=float)
    keep = (z_all >= 300.0) & (z_all <= 5.0e4) & (omega_all > 1e-8)
    z_fit = z_all[keep]
    target_omega = omega_all[keep]
    if z_fit.size < 40:
        raise RuntimeError("insufficient scalar background points for fluid matching")
    peak_weight = 1.0 + 8.0 * target_omega / np.max(target_omega)

    def residual(x: np.ndarray) -> np.ndarray:
        fluid_params = {"zc": 10.0 ** float(x[0]), "fde_zc": float(x[1])}
        pars = make_params(config, "fluid", lmax=lmax, fluid_params=fluid_params, want_cls=False)
        results = camb.get_background(pars)
        omega = np.asarray(results.get_Omega("de", z=z_fit), dtype=float)
        if np.any(~np.isfinite(omega)) or np.any(omega <= 0):
            return np.full(z_fit.size, 1e6)
        return peak_weight * (np.log(omega) - np.log(target_omega))

    fit = least_squares(
        residual,
        x0=np.array([config.log10_zc, config.fde_zc]),
        bounds=(np.array([3.2, 0.005]), np.array([4.5, 0.25])),
        xtol=1e-9,
        ftol=1e-9,
        gtol=1e-9,
        max_nfev=80,
    )
    if not fit.success:
        raise RuntimeError(f"fluid background fit failed: {fit.message}")
    params = {"zc": 10.0 ** float(fit.x[0]), "fde_zc": float(fit.x[1])}
    raw = residual(fit.x) / peak_weight
    params.update(
        {
            "log10_zc": float(fit.x[0]),
            "rms_log_omega": float(np.sqrt(np.mean(raw**2))),
            "max_abs_log_omega": float(np.max(np.abs(raw))),
            "nfev": int(fit.nfev),
        }
    )
    return params


def varied_config(config: CanonicalConfig, parameter: str, sign: int) -> CanonicalConfig:
    step = DERIVATIVE_STEPS[parameter]
    if parameter == "lnAs":
        return replace(config, As=config.As * np.exp(sign * step))
    return replace(config, **{parameter: getattr(config, parameter) + sign * step})


def generate_stock(output_dir: Path, lmax: int) -> None:
    config = CanonicalConfig()
    run_spectrum(output_dir, "lcdm", config, "lcdm", lmax)
    run_spectrum(output_dir, "scalar_n28", replace(config, n=2.8), "scalar", lmax)
    run_spectrum(output_dir, "scalar_n30", config, "scalar", lmax)
    run_spectrum(output_dir, "scalar_n32", replace(config, n=3.2), "scalar", lmax)
    for parameter in DERIVATIVE_STEPS:
        for sign, suffix in ((-1, "minus"), (1, "plus")):
            run_spectrum(
                output_dir,
                f"deriv_{parameter}_{suffix}",
                varied_config(config, parameter, sign),
                "scalar",
                lmax,
            )
    fluid = fit_effective_fluid(config, output_dir / "scalar_n30.npz", lmax)
    (output_dir / "fluid_match.json").write_text(json.dumps(json_safe(fluid), indent=2, sort_keys=True), encoding="utf-8")
    run_spectrum(output_dir, "fluid_matched", config, "fluid", lmax, fluid_params=fluid)
    (output_dir / "generation_contract.json").write_text(
        json.dumps(
            {
                "canonical_config": asdict(config),
                "zc": config.zc,
                "derivative_steps": DERIVATIVE_STEPS,
                "lmax": lmax,
                "mode": "stock",
                "nonlinear_mode": "NonLinear_none",
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def generate_nopert(output_dir: Path, lmax: int) -> None:
    config = CanonicalConfig()
    run_spectrum(output_dir, "scalar_n30_nopert", config, "scalar", lmax)
    (output_dir / "nopert_contract.json").write_text(
        json.dumps(
            {
                "canonical_config": asdict(config),
                "zc": config.zc,
                "lmax": lmax,
                "mode": "nopert",
                "nonlinear_mode": "NonLinear_none",
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("stock", "nopert"), required=True)
    parser.add_argument("--lmax", type=int, default=3000)
    args = parser.parse_args()
    if args.lmax < 300:
        raise SystemExit("lmax must be >= 300")
    if args.mode == "stock":
        generate_stock(args.output, args.lmax)
    else:
        generate_nopert(args.output, args.lmax)


if __name__ == "__main__":
    main()
