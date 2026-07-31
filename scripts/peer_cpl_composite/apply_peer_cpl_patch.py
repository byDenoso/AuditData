#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(
            f"{path}: expected exactly one match, found {count}: {old[:120]!r}"
        )
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def canonical_peer_source(path: Path) -> str:
    """Return the already-normalized composite source or fail closed.

    PeerCPL.f90 is the canonical source. The patcher must not silently rewrite it,
    because doing so made local and GitHub builds execute different code.
    """
    text = path.read_text(encoding="utf-8")
    required = (
        "module PeerCPL",
        "type, extends(TEarlyQuintessence) :: TPeerCPL",
        "call this%TQuintessence%Init(State)",
        "this%State%grho_no_de(a)",
        "if (a <= 0._dl) then",
        "end subroutine TPeerCPL_Effective_w_wa",
        "end module PeerCPL",
    )
    forbidden = (
        "end suboutine",
        "end module PeerCPLS",
        "CompositeState",
        "wtot",
        "weff",
        "waeff",
    )
    missing = [token for token in required if token not in text]
    surviving = [token for token in forbidden if token in text]
    if missing or surviving:
        raise RuntimeError(
            f"{path}: non-canonical PeerCPL source; missing={missing}, forbidden={surviving}"
        )
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

    # Copy the exact reviewed source. No hidden source normalization is allowed.
    (fortran / "PeerCPL.f90").write_text(
        canonical_peer_source(peer_source), encoding="utf-8"
    )

    replace_once(
        fortran / "Makefile_main",
        "DARKENERGY_FILES  ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize DarkEnergyQuintessence",
        "DARKENERGY_FILES  ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize DarkEnergyQuintessence PeerCPL",
    )

    patch_dynamic_quintessence_dispatch(quintessence)
    # The derived composite EvolveBackground needs access to CAMB's state pointer.
    replace_once(
        quintessence,
        "class(CAMBdata), pointer, private :: State",
        "class(CAMBdata), pointer :: State",
    )

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
    if "__composite_state" in fragment:
        raise RuntimeError("Python fragment still exposes removed CompositeState pointer")
    text = text.replace(marker, fragment + marker, 1)
    old_map = 'F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF})'
    new_map = 'F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF, "peer_cpl": PeerCPL})'
    if text.count(old_map) != 1:
        raise RuntimeError("dark_energy.py class-name map mismatch")
    py.write_text(text.replace(old_map, new_map, 1), encoding="utf-8")

    patched_peer = (fortran / "PeerCPL.f90").read_text(encoding="utf-8")
    if patched_peer != canonical_peer_source(peer_source):
        raise RuntimeError("copied PeerCPL source differs from canonical source")

    patched_quintessence = quintessence.read_text(encoding="utf-8")
    if patched_quintessence.count("EvolveBackgroundDispatch") != 4:
        raise RuntimeError("dynamic quintessence dispatch patch did not apply exactly")
    if "class(CAMBdata), pointer, private :: State" in patched_quintessence:
        raise RuntimeError("Quintessence State pointer remains private")

    print(f"patched CAMB source at {root}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    args = parser.parse_args()
    apply(args.source_root)


if __name__ == "__main__":
    main()
