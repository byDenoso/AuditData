#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected exactly one match, found {count}: {old[:120]!r}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def replace_source_once(text: str, old: str, new: str, *, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"PeerCPL source: expected one {label}, found {count}")
    return text.replace(old, new, 1)


def normalized_peer_source(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    repairs = {
        "end suboutine TPeerCPL_Effective_w_wa": "end subroutine TPeerCPL_Effective_w_wa",
        "end module PeerCPLS": "end module PeerCPL",
        "this%TQuintessence%BackgroundDensityAndPressure": (
            "this%TEarlyQuintessence%BackgroundDensityAndPressure"
        ),
        "this%TQuintessence%PerturbationEvolve": (
            "this%TEarlyQuintessence%PerturbationEvolve"
        ),
    }
    for old, new in repairs.items():
        text = replace_source_once(text, old, new, label=old)

    text = replace_source_once(
        text,
        "subroutine TPeerCPL_EvolveBackground(this, num, a, y, yprime)\n"
        "    ! Exact scalar-background evolution in the presence of the late CPL density.\n"
        "    ! Variables are phi=y(1), a^2 phi'=y(2), matching TQuintessence.\n"
        "    class(TPeerCPL), intent(in) :: this\n"
        "    integer, intent(in) :: num\n"
        "    real(dl), intent(in) :: a, y(num)\n"
        "    real(dl), intent(out) :: yprime(num)",
        "subroutine TPeerCPL_EvolveBackground(this, num, a, y, yprime)\n"
        "    ! Exact scalar-background evolution in the presence of the late CPL density.\n"
        "    ! Variables are phi=y(1), a^2 phi'=y(2), matching TQuintessence.\n"
        "    class(TPeerCPL) :: this\n"
        "    integer :: num\n"
        "    real(dl) :: a, y(num), yprime(num)",
        label="EvolveBackground override interface",
    )
    text = replace_source_once(
        text,
        "subroutine PeerPerturbations(this, a, k, y, w_ix, dgrho_peer, dgq_peer)\n"
        "    class(TPeerCPL), intent(in) :: this\n"
        "    real(dl), intent(in) :: a, k\n"
        "    real(dl), intent(in) :: y(:)",
        "subroutine PeerPerturbations(this, a, k, y, w_ix, dgrho_peer, dgq_peer)\n"
        "    class(TPeerCPL), intent(in) :: this\n"
        "    real(dl), intent(in) :: a, k\n"
        "    real(dl), intent(in) :: y(*)",
        label="assumed-size perturbation helper",
    )

    forbidden = (
        "end suboutine",
        "end module PeerCPLS",
        "this%TQuintessence%BackgroundDensityAndPressure",
        "this%TQuintessence%PerturbationEvolve",
    )
    if any(token in text for token in forbidden):
        raise RuntimeError(f"{path}: invalid source token survived normalization")
    return text


def patch_dynamic_quintessence_dispatch(path: Path) -> None:
    wrapper_anchor = (
        "    end subroutine EvolveBackground\n\n\n"
        "    real(dl) function TQuintessence_phidot_start(this,phi)"
    )
    wrapper = (
        "    end subroutine EvolveBackground\n\n"
        "    subroutine EvolveBackgroundDispatch(this,num,a,y,yprime)\n"
        "    ! Preserve dynamic dispatch for derived quintessence models during\n"
        "    ! the linear-a integration phase of TEarlyQuintessence_Init.\n"
        "    class(TQuintessence) :: this\n"
        "    integer num\n"
        "    real(dl) y(num), yprime(num), a\n\n"
        "    call this%EvolveBackground(num, a, y, yprime)\n"
        "    end subroutine EvolveBackgroundDispatch\n\n\n"
        "    real(dl) function TQuintessence_phidot_start(this,phi)"
    )
    replace_once(path, wrapper_anchor, wrapper)
    replace_once(
        path,
        "        call dverk(this,NumEqs,EvolveBackground,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)\n"
        "        if (.not. this%check_error(afrom, aend)) return\n"
        "        call EvolveBackground(this,NumEqs,aend,y,w(:,1))",
        "        call dverk(this,NumEqs,EvolveBackgroundDispatch,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)\n"
        "        if (.not. this%check_error(afrom, aend)) return\n"
        "        call EvolveBackgroundDispatch(this,NumEqs,aend,y,w(:,1))",
    )


def apply(root: Path) -> None:
    root = root.resolve()
    here = Path(__file__).resolve().parent
    fortran = root / "fortran"
    py = root / "camb" / "dark_energy.py"
    peer_source = here / "PeerCPL.f90"
    quintessence = fortran / "DarkEnergyQuintessence.f90"
    required = [
        fortran / "DarkEnergyInterface.f90",
        fortran / "DarkEnergyPPF.f90",
        quintessence,
        fortran / "equations.f90",
        fortran / "Makefile_main",
        py,
        peer_source,
        here / "peer_cpl_python.pyfrag",
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    (fortran / "PeerCPL.f90").write_text(
        normalized_peer_source(peer_source), encoding="utf-8"
    )

    replace_once(
        fortran / "Makefile_main",
        "DARKENERGY_FILES  ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize DarkEnergyQuintessence",
        "DARKENERGY_FILES  ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize DarkEnergyQuintessence PeerCPL",
    )

    patch_dynamic_quintessence_dispatch(quintessence)

    replace_once(
        fortran / "DarkEnergyInterface.f90",
        "function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &",
        "function diff_rhopi_Add_Term(this, a, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &",
    )
    replace_once(
        fortran / "DarkEnergyInterface.f90",
        "real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &",
        "real(dl), intent(in) :: a, dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &",
    )

    replace_once(
        fortran / "DarkEnergyPPF.f90",
        "function TDarkEnergyPPF_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &",
        "function TDarkEnergyPPF_diff_rhopi_Add_Term(this, a, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &",
    )
    replace_once(
        fortran / "DarkEnergyPPF.f90",
        "real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &",
        "real(dl), intent(in) :: a, dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &",
    )

    replace_once(
        fortran / "equations.f90",
        "State%CP%DarkEnergy%diff_rhopi_Add_Term(dgrho_de, dgq_de, grho, &",
        "State%CP%DarkEnergy%diff_rhopi_Add_Term(a, dgrho_de, dgq_de, grho, &",
    )

    text = py.read_text(encoding="utf-8")
    marker = "\n\n# short names for models that support w/wa\n"
    if marker not in text:
        raise RuntimeError("dark_energy.py insertion marker missing")
    fragment = (here / "peer_cpl_python.pyfrag").read_text(encoding="utf-8")
    text = text.replace(marker, fragment + marker, 1)
    old_map = 'F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF})'
    new_map = 'F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF, "peer_cpl": PeerCPL})'
    if text.count(old_map) != 1:
        raise RuntimeError("dark_energy.py class-name map mismatch")
    py.write_text(text.replace(old_map, new_map, 1), encoding="utf-8")

    patched_peer = (fortran / "PeerCPL.f90").read_text(encoding="utf-8")
    required_tokens = (
        "module PeerCPL",
        "type, extends(TEarlyQuintessence) :: TPeerCPL",
        "this%TEarlyQuintessence%BackgroundDensityAndPressure",
        "this%TEarlyQuintessence%PerturbationEvolve",
        "end subroutine TPeerCPL_Effective_w_wa",
        "end module PeerCPL",
    )
    for token in required_tokens:
        if token not in patched_peer:
            raise RuntimeError(f"patched PeerCPL source missing {token!r}")

    patched_quintessence = quintessence.read_text(encoding="utf-8")
    if patched_quintessence.count("EvolveBackgroundDispatch") != 4:
        raise RuntimeError("dynamic quintessence dispatch patch did not apply exactly")

    print(f"patched CAMB source at {root}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    args = parser.parse_args()
    apply(args.source_root)


if __name__ == "__main__":
    main()
