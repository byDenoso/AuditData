from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PATCH_PATH = ROOT / "scripts" / "peer_cpl_composite" / "apply_peer_cpl_patch.py"
SPEC = importlib.util.spec_from_file_location("apply_peer_cpl_patch", PATCH_PATH)
assert SPEC and SPEC.loader
PATCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATCH)


class PeerCPLPatchTests(unittest.TestCase):
    def make_source_tree(self, root: Path, *, drift_makefile: bool = False) -> Path:
        fortran = root / "fortran"
        camb = root / "camb"
        fortran.mkdir(parents=True)
        camb.mkdir(parents=True)

        make_line = (
            "DARKENERGY_FILES  ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize "
            "DarkEnergyQuintessence"
        )
        if drift_makefile:
            make_line += " UnexpectedModel"
        (fortran / "Makefile_main").write_text(make_line + "\n", encoding="utf-8")

        (fortran / "DarkEnergyInterface.f90").write_text(
            "function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &\n"
            "real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &\n",
            encoding="utf-8",
        )
        (fortran / "DarkEnergyPPF.f90").write_text(
            "function TDarkEnergyPPF_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &\n"
            "real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &\n",
            encoding="utf-8",
        )
        (fortran / "equations.f90").write_text(
            "State%CP%DarkEnergy%diff_rhopi_Add_Term(dgrho_de, dgq_de, grho, &\n",
            encoding="utf-8",
        )
        (fortran / "DarkEnergyQuintessence.f90").write_text(
            "    end subroutine EvolveBackground\n\n\n"
            "    real(dl) function TQuintessence_phidot_start(this,phi)\n"
            "        call dverk(this,NumEqs,EvolveBackground,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)\n"
            "        if (.not. this%check_error(afrom, aend)) return\n"
            "        call EvolveBackground(this,NumEqs,aend,y,w(:,1))\n",
            encoding="utf-8",
        )
        (camb / "dark_energy.py").write_text(
            "class EarlyQuintessence:\n    pass\n"
            "\n\n# short names for models that support w/wa\n"
            'F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF})\n',
            encoding="utf-8",
        )
        return root

    def test_patch_applies_all_composite_hooks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = self.make_source_tree(Path(tmp))
            PATCH.apply(source)

            peer = (source / "fortran" / "PeerCPL.f90").read_text(encoding="utf-8")
            quint = (source / "fortran" / "DarkEnergyQuintessence.f90").read_text(
                encoding="utf-8"
            )
            interface = (source / "fortran" / "DarkEnergyInterface.f90").read_text(
                encoding="utf-8"
            )
            equations = (source / "fortran" / "equations.f90").read_text(encoding="utf-8")
            py = (source / "camb" / "dark_energy.py").read_text(encoding="utf-8")

            self.assertIn("end subroutine TPeerCPL_Effective_w_wa", peer)
            self.assertIn("end module PeerCPL", peer)
            self.assertNotIn("end suboutine", peer)
            self.assertIn("this%TEarlyQuintessence%BackgroundDensityAndPressure", peer)
            self.assertIn("this%TEarlyQuintessence%PerturbationEvolve", peer)
            self.assertIn("real(dl), intent(in) :: y(*)", peer)
            self.assertEqual(quint.count("EvolveBackgroundDispatch"), 4)
            self.assertIn("diff_rhopi_Add_Term(this, a,", interface)
            self.assertIn("diff_rhopi_Add_Term(a, dgrho_de", equations)
            self.assertIn("class PeerCPL", py)
            self.assertIn('"peer_cpl": PeerCPL', py)

    def test_patch_fails_closed_on_upstream_source_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = self.make_source_tree(Path(tmp), drift_makefile=True)
            with self.assertRaises(RuntimeError):
                PATCH.apply(source)


if __name__ == "__main__":
    unittest.main()
