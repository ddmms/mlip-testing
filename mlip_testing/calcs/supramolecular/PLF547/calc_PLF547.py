"""Run calculations for PLF547 protein-ligand benchmark."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from ase import Atoms
from ase.io import write
import mlipx
from mlipx.abc import NodeWithCalculator
import numpy as np
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


class PLF547Benchmark(zntrack.Node):
    """
    Benchmark model for PLF547 protein-ligand interactions.

    Evaluates interaction energies for protein fragment-ligand interactions.
    Each system consists of complex, protein, and ligand fragments.
    Computes interaction energy = E(complex) - E(protein) - E(ligand)
    """

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def extract_charge_and_selections(
        pdb_path: Path,
    ) -> tuple[float, float, float, str, str]:
        """
        Extract charge and selection information from PDB REMARK lines.

        Parameters
        ----------
        pdb_path : Path
            Path to the PDB file.

        Returns
        -------
        tuple[float, float, float, str, str]
            Total charge, charge A, charge B, selection A, selection B.
        """
        total_charge = qa = qb = 0.0
        selection_a = selection_b = ""

        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("REMARK"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        break
                    continue

                parts = line.split()
                if len(parts) < 3:
                    continue

                tag = parts[1].lower()

                if tag == "charge":
                    total_charge = float(parts[2])
                elif tag == "charge_a":
                    qa = float(parts[2])
                elif tag == "charge_b":
                    qb = float(parts[2])
                elif tag == "selection_a":
                    selection_a = " ".join(parts[2:])
                elif tag == "selection_b":
                    selection_b = " ".join(parts[2:])

        return total_charge, qa, qb, selection_a, selection_b

    @staticmethod
    def separate_protein_ligand_simple(pdb_path: Path):
        """
        Separate protein and ligand based on residue names.

        Parameters
        ----------
        pdb_path : Path
            Path to the PDB file.

        Returns
        -------
        tuple
            All atoms, protein atoms, ligand atoms.
        """
        import MDAnalysis as MDAnalysis

        # Load with MDAnalysis
        u = MDAnalysis.Universe(str(pdb_path))

        # Simple separation: ligand = UNK residues, protein = everything else
        protein_atoms = []
        ligand_atoms = []

        for atom in u.atoms:
            if atom.resname.strip().upper() in ["UNK", "LIG", "MOL"]:
                ligand_atoms.append(atom)
            else:
                protein_atoms.append(atom)

        return u.atoms, protein_atoms, ligand_atoms

    @staticmethod
    def mda_atoms_to_ase(atom_list, charge: float, identifier: str) -> Atoms:
        """
        Convert MDAnalysis atoms to ASE Atoms object.

        Parameters
        ----------
        atom_list : list
            List of MDAnalysis atoms.
        charge : float
            Charge of the system.
        identifier : str
            System identifier.

        Returns
        -------
        Atoms
            ASE Atoms object.
        """
        if not atom_list:
            atoms = Atoms()
            atoms.info.update({"charge": charge, "identifier": identifier})
            return atoms

        symbols = []
        positions = []

        for atom in atom_list:
            # Get element symbol
            try:
                elem = (atom.element or "").strip().title()
            except AttributeError:
                elem = ""

            if not elem:
                # Fallback: first letter of atom name
                elem = "".join([c for c in atom.name if c.isalpha()])[:1].title() or "C"

            symbols.append(elem)
            positions.append(atom.position)

        atoms = Atoms(symbols=symbols, positions=np.array(positions))
        atoms.info.update({"charge": int(round(charge)), "identifier": identifier})
        return atoms

    @staticmethod
    def process_pdb_file(pdb_path: Path) -> dict[str, Atoms]:
        """
        Generate complex and separated fragments from one PDB file.

        Parameters
        ----------
        pdb_path : Path
            Path to the PDB file.

        Returns
        -------
        dict[str, Atoms]
            Dictionary with 'complex', 'protein', and 'ligand' fragments.
        """
        total_charge, charge_a, charge_b, _, _ = (
            PLF547Benchmark.extract_charge_and_selections(pdb_path)
        )

        try:
            all_atoms, protein_atoms, ligand_atoms = (
                PLF547Benchmark.separate_protein_ligand_simple(pdb_path)
            )

            if len(ligand_atoms) == 0:
                logging.warning(f"No ligand atoms found in {pdb_path.name}")
                return {}

            if len(protein_atoms) == 0:
                logging.warning(f"No protein atoms found in {pdb_path.name}")
                return {}

            base_id = pdb_path.stem

            complex_atoms = PLF547Benchmark.mda_atoms_to_ase(
                list(all_atoms), total_charge, base_id
            )
            protein_frag = PLF547Benchmark.mda_atoms_to_ase(
                protein_atoms, charge_a, base_id
            )
            ligand = PLF547Benchmark.mda_atoms_to_ase(ligand_atoms, charge_b, base_id)

            return {"complex": complex_atoms, "protein": protein_frag, "ligand": ligand}

        except Exception as e:
            logging.warning(f"Error processing {pdb_path}: {e}")
            return {}

    @staticmethod
    def parse_plf547_references(path: Path) -> dict[str, float]:
        """
        Parse PLF547 reference interaction energies (kcal/mol -> eV).

        Parameters
        ----------
        path : Path
            Path to the reference energies file.

        Returns
        -------
        dict[str, float]
            Dictionary mapping system names to reference energies in eV.
        """
        ref: dict[str, float] = {}

        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.lower().startswith("no.") or line.startswith("-"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            try:
                energy_kcal = float(parts[-1])
            except ValueError:
                continue

            # Extract full identifier with residue type
            full_identifier = parts[1].replace(".pdb", "")

            # Extract base identifier by removing residue type suffix
            # Format: "2P4Y_28_met" -> "2P4Y_28"
            identifier_parts = full_identifier.split("_")
            if len(identifier_parts) >= 3:
                # Assume last part is residue type (met, arg, bbn, etc.)
                base_identifier = "_".join(identifier_parts[:-1])
            else:
                # Fallback: use full identifier if format is unexpected
                base_identifier = full_identifier

            energy_ev = energy_kcal * KCAL_TO_EV  # Convert to eV
            ref[base_identifier] = energy_ev

        return ref

    def run(self):
        """Run PLF547 benchmark calculations."""
        calc = self.model.get_calculator()

        # Get benchmark data
        base_dir = (
            get_benchmark_data("protein-ligand-data_PLA15_PLF547.zip")
            / "protein-ligand-data_PLA15_PLF547"
        )

        # PLF547 specific paths
        plf547_dir = base_dir / "PLF547_pdbs"
        plf547_ref_file = plf547_dir / "reference_energies.txt"

        if not plf547_dir.exists() or not plf547_ref_file.exists():
            logging.warning("PLF547 data directory or reference file not found")
            return

        # Write output structures
        write_dir = OUT_PATH / self.model_name
        write_dir.mkdir(parents=True, exist_ok=True)

        plf547_refs = self.parse_plf547_references(plf547_ref_file)
        pdb_files = list(plf547_dir.glob("*.pdb"))

        results = []

        for pdb_file in tqdm(pdb_files, desc="PLF547"):
            identifier = pdb_file.stem
            if identifier not in plf547_refs:
                continue
            fragments = self.process_pdb_file(pdb_file)
            if not fragments:
                continue
            try:
                # Direct energy computation
                fragments["complex"].calc = calc
                e_complex = fragments["complex"].get_potential_energy()
                fragments["protein"].calc = calc
                e_protein = fragments["protein"].get_potential_energy()
                fragments["ligand"].calc = calc
                e_ligand = fragments["ligand"].get_potential_energy()
                e_int_model = e_complex - e_protein - e_ligand
                e_int_ref = plf547_refs[identifier]

                results.append(
                    {
                        "identifier": identifier,
                        "E_int_ref": e_int_ref,
                        f"E_int_{self.model_name}": e_int_model,
                        "error_eV": e_int_model - e_int_ref,
                        "error_kcal": (e_int_model - e_int_ref) * EV_TO_KCAL,
                        "complex_atoms": len(fragments["complex"]),
                        "protein_atoms": len(fragments["protein"]),
                        "ligand_atoms": len(fragments["ligand"]),
                    }
                )

                # Store additional info in atoms and save individual structure
                fragments["complex"].info.update(
                    {
                        "identifier": identifier,
                        "model": self.model_name,
                        f"E_int_{self.model_name}": e_int_model,
                        "E_int_ref": e_int_ref,
                        "error_eV": e_int_model - e_int_ref,
                        "error_kcal": (e_int_model - e_int_ref) * EV_TO_KCAL,
                    }
                )

                # Save individual structure file with identifier as filename
                write(
                    write_dir / f"{identifier}.xyz",
                    fragments["complex"],
                    format="extxyz",
                )

            except Exception as e:
                logging.warning(f"Error processing PLF547 {identifier}: {e}")
                continue

        # Save results
        if results:
            plf547_df = pd.DataFrame(results)
            plf547_df.to_csv(write_dir / "plf547_results.csv", index=False)
            plf547_mae = float(np.abs(plf547_df["error_kcal"]).mean())

            mae_data = {"PLF547_MAE_kcal": plf547_mae}
            print(f"PLF547 MAE: {plf547_mae}")
            with open(write_dir / "mae_results.json", "w") as f:
                json.dump(mae_data, f, indent=2)


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
            benchmark = PLF547Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_plf547():
    """Run PLF547 benchmark via pytest."""
    build_project(repro=True)


if __name__ == "__main__":
    build_project()
