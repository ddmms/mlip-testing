"""Run calculations for PLF547 benchmark."""

from __future__ import annotations

from copy import deepcopy
import logging
from pathlib import Path

from ase import Atoms, units
from ase.calculators.calculator import Calculator
from ase.io import write
import mlipx
from mlipx.abc import NodeWithCalculator
import numpy as np
from tqdm import tqdm
import zntrack

from mlip_testing.calcs.models.models import MODELS
from mlip_testing.calcs.utils.utils import chdir, get_benchmark_data

# Local directory to store output data
OUT_PATH = Path(__file__).parent / "outputs"

# Constants
KCAL_PER_MOL_TO_EV = units.kcal / units.mol
EV_TO_KCAL_PER_MOL = 1.0 / KCAL_PER_MOL_TO_EV


class PLF547Benchmark(zntrack.Node):
    """
    Benchmark model for PLF547 dataset.

    Protein-ligand fragment interaction energy calculations.
    - 547 protein fragment-ligand interactions
    - Each system consists of complex, protein fragment, and ligand
    - Computes interaction energy = E(complex) - E(protein) - E(ligand)
    """

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def compute_interaction_energy(
        complex_e: float, protein_e: float, ligand_e: float
    ) -> float:
        """
        Compute interaction energy.

        Parameters
        ----------
        complex_e
            Energy of the complex.
        protein_e
            Energy of the protein fragment.
        ligand_e
            Energy of the ligand.

        Returns
        -------
        float
            Interaction energy.
        """
        return complex_e - protein_e - ligand_e

    @staticmethod
    def extract_charge_and_selections(
        pdb_path: Path,
    ) -> tuple[float, float, float, str, str]:
        """
        Extract charge and selection information from PDB REMARK lines.

        Parameters
        ----------
        pdb_path
            Path to PDB file to parse.

        Returns
        -------
        tuple[float, float, float, str, str]
            Total charge, charge_a, charge_b, selection_a, selection_b.
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
        Separate protein and ligand atoms based on residue names.

        Parameters
        ----------
        pdb_path
            Path to PDB file to process.

        Returns
        -------
        tuple
            All atoms, protein atoms, ligand atoms from MDAnalysis.
        """
        import MDAnalysis

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
        atom_list
            List of MDAnalysis atoms to convert.
        charge
            Total charge of the system.
        identifier
            System identifier string.

        Returns
        -------
        Atoms
            ASE Atoms object with charge and identifier in info dict.
        """
        if not atom_list:
            atoms = Atoms()
            atoms.info.update({"charge": int(round(charge)), "identifier": identifier})
            return atoms

        symbols = []
        positions = []

        for atom in atom_list:
            # Get element symbol
            try:
                elem = (atom.element or "").strip().title()
            except (AttributeError, TypeError):
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
        Parse one PDB file and return complex and separated fragments.

        Parameters
        ----------
        pdb_path
            Path to PDB file to process.

        Returns
        -------
        dict[str, Atoms]
            Dictionary containing 'complex', 'protein', and 'ligand' ASE Atoms objects.
        """
        (
            total_charge,
            charge_a,
            charge_b,
            _,
            _,
        ) = PLF547Benchmark.extract_charge_and_selections(pdb_path)

        try:
            (
                all_atoms,
                protein_atoms,
                ligand_atoms,
            ) = PLF547Benchmark.separate_protein_ligand_simple(pdb_path)

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

        except (OSError, ValueError, ImportError) as e:
            logging.warning(f"Error processing {pdb_path}: {e}")
            return {}

    @staticmethod
    def parse_plf547_references(path: Path) -> dict[str, float]:
        """
        Parse PLF547 reference interaction energies (kcal/mol -> eV).

        Parameters
        ----------
        path
            Path to reference energies file.

        Returns
        -------
        dict[str, float]
            Dictionary mapping system identifiers to reference energies in eV.
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

            # Extract base identifier by removing residue type suffix
            # Format: "2P4Y_28_met" -> "2P4Y_28"
            full_identifier = parts[1].replace(".pdb", "")
            identifier_parts = full_identifier.split("_")
            if len(identifier_parts) >= 3:
                # Assume last part is residue type (met, arg, bbn, etc.)
                base_identifier = "_".join(identifier_parts[:-1])
            else:
                # Fallback: use full identifier if format is unexpected
                base_identifier = full_identifier

            energy_ev = energy_kcal * KCAL_PER_MOL_TO_EV  # Convert to eV
            ref[base_identifier] = energy_ev

        return ref

    @staticmethod
    def evaluate_energies(fragments: dict[str, Atoms], calc: Calculator) -> None:
        """
        Evaluate energies for PLF547 structures.

        Parameters
        ----------
        fragments
            Dictionary containing complex, protein, and ligand atoms.
        calc
            Calculator to use to evaluate structure energy.
        """
        for atoms in fragments.values():
            atoms.calc = deepcopy(calc)
            atoms.get_potential_energy()

    def run(self):
        """Run PLF547 energy calculations."""
        calc = self.model.get_calculator()
        base_dir = (
            get_benchmark_data("protein-ligand-data_PLA15_PLF547.zip")
            / "protein-ligand-data_PLA15_PLF547"
        )
        plf547_dir = base_dir / "PLF547_pdbs"
        plf547_ref_file = plf547_dir / "reference_energies.txt"

        plf547_refs = self.parse_plf547_references(plf547_ref_file)
        pdb_files = list(plf547_dir.glob("*.pdb"))

        complex_atoms_list = []

        for pdb_file in tqdm(
            pdb_files, desc=f"Processing PLF547 for model: {self.model_name}"
        ):
            identifier = pdb_file.stem
            if identifier not in plf547_refs:
                continue

            fragments = self.process_pdb_file(pdb_file)
            if not fragments:
                continue

            try:
                # Evaluate with the model
                self.evaluate_energies(fragments, calc)

                # Calculate interaction energies
                complex_e = fragments["complex"].get_potential_energy()
                protein_e = fragments["protein"].get_potential_energy()
                ligand_e = fragments["ligand"].get_potential_energy()
                pred_int_energy = self.compute_interaction_energy(
                    complex_e, protein_e, ligand_e
                )

                ref_int_energy = plf547_refs[identifier]

                # Calculate errors
                error_ev = pred_int_energy - ref_int_energy
                error_kcal = error_ev * EV_TO_KCAL_PER_MOL

                # Store information in complex atoms
                complex_atoms = fragments["complex"]
                complex_atoms.info["model"] = self.model_name
                complex_atoms.info["interaction_energy"] = pred_int_energy
                complex_atoms.info["ref_interaction_energy"] = ref_int_energy
                complex_atoms.info["error_ev"] = error_ev
                complex_atoms.info["error_kcal"] = error_kcal
                complex_atoms.info["sys_id"] = identifier
                complex_atoms.info["system"] = identifier

                complex_atoms_list.append(complex_atoms)

            except (RuntimeError, ValueError, KeyError) as e:
                logging.warning(f"Error calculating energies for {identifier}: {e}")
                continue

        # Write output structures
        write_dir = OUT_PATH / self.model_name
        write_dir.mkdir(parents=True, exist_ok=True)

        # Save individual complex atoms files for each system
        for i, atoms in enumerate(complex_atoms_list):
            # Clear calculator to avoid array broadcasting issues
            atoms_copy = atoms.copy()
            atoms_copy.calc = None

            # Write each system to its own file
            system_file = write_dir / f"{i}.xyz"
            write(system_file, atoms_copy, format="extxyz")


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
