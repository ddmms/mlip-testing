"""Run calculations for S30L benchmark."""

from __future__ import annotations

import json
from pathlib import Path
import warnings

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.io import read, write
import mlipx
from mlipx.abc import NodeWithCalculator
import pandas as pd
from tqdm import tqdm
import zntrack

from mlip_testing.calcs.models.models import MODELS
from mlip_testing.calcs.utils.utils import chdir, get_benchmark_data

# Local directory to store output data
OUT_PATH = Path(__file__).parent / "outputs"

# Constants
KCAL_TO_EV = 0.04336414
EV_TO_KCAL = 1.0 / KCAL_TO_EV


class S30LBenchmark(zntrack.Node):
    """
    Benchmark model for S30L dataset.

    Evaluates interaction energies for 30 large host-guest complexes.
    Each system consists of host (A), guest (B), and complex (AB).

    Computes interaction energy = E(complex) - E(host) - E(guest)
    """

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def _read_charge(folder: Path) -> float:
        """
        Read charge from .CHRG file in folder.

        Parameters
        ----------
        folder : Path
            Path to the folder containing .CHRG file.

        Returns
        -------
        float
            Charge value, 0.0 if not found.
        """
        for f in folder.iterdir():
            if f.name.upper() == ".CHRG":
                try:
                    return float(f.read_text().strip())
                except ValueError:
                    warnings.warn(
                        f"Invalid charge in {f} - assuming neutral.", stacklevel=2
                    )
        return 0.0

    @staticmethod
    def _read_atoms(folder: Path, ident: str) -> Atoms:
        """
        Read atoms from coord file in folder.

        Parameters
        ----------
        folder : Path
            Path to the folder containing coord file.
        ident : str
            System identifier.

        Returns
        -------
        Atoms
            ASE Atoms object.
        """
        coord = next(
            (p for p in folder.iterdir() if p.name.lower().startswith("coord")), None
        )
        if coord is None:
            raise FileNotFoundError(f"No coord file in {folder}")
        atoms = read(coord, format="turbomole")
        atoms.info.update(
            {"identifier": ident, "charge": int(S30LBenchmark._read_charge(folder))}
        )
        return atoms

    @staticmethod
    def load_complex(index: int, root: Path) -> dict[str, Atoms]:
        """
        Load host, guest, and complex structures for given index.

        Parameters
        ----------
        index : int
            System index number.
        root : Path
            Root directory containing system data.

        Returns
        -------
        dict[str, Atoms]
            Dictionary with 'host', 'guest', and 'complex' Atoms objects.
        """
        base = root / f"{index}"
        if not base.exists():
            raise FileNotFoundError(base)
        return {
            "host": S30LBenchmark._read_atoms(base / "A", f"{index}_host"),
            "guest": S30LBenchmark._read_atoms(base / "B", f"{index}_guest"),
            "complex": S30LBenchmark._read_atoms(base / "AB", f"{index}_complex"),
        }

    @staticmethod
    def interaction_energy(frags: dict[str, Atoms], calc: Calculator) -> float:
        """
        Calculate interaction energy from fragments.

        Parameters
        ----------
        frags : dict[str, Atoms]
            Dictionary containing 'complex', 'host', and 'guest' fragments.
        calc : Calculator
            ASE calculator for energy calculations.

        Returns
        -------
        float
            Interaction energy in eV.
        """
        frags["complex"].calc = calc
        e_complex = frags["complex"].get_potential_energy()
        frags["host"].calc = calc
        e_host = frags["host"].get_potential_energy()
        frags["guest"].calc = calc
        e_guest = frags["guest"].get_potential_energy()
        return e_complex - e_host - e_guest

    @staticmethod
    def parse_references(path: Path) -> dict[int, float]:
        """
        Parse reference energies from text file.

        Parameters
        ----------
        path : Path
            Path to the reference energies file.

        Returns
        -------
        dict[int, float]
            Dictionary mapping system indices to reference energies in eV.
        """
        refs: dict[int, float] = {}
        for idx, ln in enumerate(path.read_text().splitlines()):
            ln = ln.strip()
            if not ln:
                continue
            kcal = float(ln.split()[0])
            refs[idx + 1] = kcal * KCAL_TO_EV
        return refs

    def run(self):
        """Run S30L benchmark calculations."""
        calc = self.model.get_calculator()

        # Get benchmark data
        base_dir = get_benchmark_data("S30L.zip") / "S30L/s30l_test_set"
        ref_file = base_dir / "references_s30.txt"
        refs = self.parse_references(ref_file)

        rows = []
        complex_atoms_list = []

        for idx in tqdm(range(1, 31), desc=f"Benchmarking {self.model_name}"):
            try:
                fragments = self.load_complex(idx, base_dir)
                e_model = self.interaction_energy(fragments, calc)
                e_ref = refs[idx]
                error_ev = e_model - e_ref
                error_kcal = error_ev * EV_TO_KCAL

                rows.append(
                    {
                        "Index": idx,
                        "E_ref (eV)": e_ref,
                        f"E_{self.model_name} (eV)": e_model,
                        "Error (eV)": error_ev,
                        "Error (kcal/mol)": error_kcal,
                        "n_atoms": len(fragments["complex"]),
                    }
                )

                # Store additional info in complex atoms
                fragments["complex"].info.update(
                    {
                        "Index": idx,
                        "model": self.model_name,
                        "E_model_eV": e_model,
                        "E_ref_eV": e_ref,
                        "error_eV": error_ev,
                        "error_kcal": error_kcal,
                    }
                )
                complex_atoms_list.append(fragments["complex"])

            except Exception as e:
                print(f"Error processing system {idx}: {e}")
                continue

        # Create results DataFrame
        df = pd.DataFrame(rows)

        # Write output structures and results
        write_dir = OUT_PATH / self.model_name
        write_dir.mkdir(parents=True, exist_ok=True)

        # Write individual structure files for each system
        if complex_atoms_list:
            for i, atoms in enumerate(complex_atoms_list):
                # Create fresh atoms object to avoid array contamination issues
                fresh_atoms = Atoms(
                    symbols=atoms.get_chemical_symbols(),
                    positions=atoms.positions.copy(),
                    cell=atoms.cell.copy() if atoms.cell is not None else None,
                    pbc=atoms.pbc.copy() if atoms.pbc is not None else False,
                )
                fresh_atoms.info.update(atoms.info)

                # Write each system to its own file
                system_file = write_dir / f"{i}.xyz"
                write(system_file, fresh_atoms, format="extxyz")

        if not df.empty:
            df.to_csv(write_dir / "s30l_results.csv", index=False)

            # Calculate and save MAE
            mae = df["Error (kcal/mol)"].abs().mean()
            mae_data = {"MAE_kcal": float(mae)}

            with open(write_dir / "mae_results.json", "w") as f:
                json.dump(mae_data, f, indent=2)

            print(f"MAE for {self.model_name}: {mae:.2f} kcal/mol")


def build_project(repro: bool = False) -> None:
    """
    Build mlipx project.

    Parameters
    ----------
    repro
        Whether to call dvc repro -f after building.
    """
    project = mlipx.Project()
    benchmark_node_dict = {}

    for model_name, model in MODELS.items():
        with project.group(model_name):
            benchmark = S30LBenchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_s30l():
    """Run S30L benchmark via pytest."""
    build_project(repro=True)


if __name__ == "__main__":
    build_project()
