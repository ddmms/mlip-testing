"""Run calculations for OC157 benchmark."""

from __future__ import annotations

from copy import deepcopy
import pathlib
from pathlib import Path
import zipfile

from ase.io import read, write
import mlipx
from mlipx.abc import NodeWithCalculator
import requests
from tqdm import tqdm
import zntrack

from mlip_testing.calcs.models.models import MODELS

# Local directory to store output data
OUT_PATH = Path(__file__).parent / "outputs"

# GitHub API URL (for listing contents if needed)
BENCHMARK_DATA_URL = (
    "https://api.github.com/repos/joehart2001/mlipx/contents/benchmark_data"
)

# Raw download URL (for direct downloading)
BENCHMARK_DATA_DOWNLOAD_URL = (
    "https://raw.githubusercontent.com/joehart2001/mlipx/main/benchmark_data/"
)

# Local cache directory
BENCHMARK_DATA_DIR = pathlib.Path.home() / ".cache" / "my_benchmark"


def get_benchmark_data(name: str, force: bool = False) -> Path:
    """
    Retrieve benchmark data.

    If it's a .zip, download and extract it.

    Parameters
    ----------
    name
        Name of benchmark data file.
    force
        Whether to ignore cached download.

    Returns
    -------
    Path
        Path to extracted data.
    """
    uri = f"{BENCHMARK_DATA_DOWNLOAD_URL}/{name}"
    local_path = Path(BENCHMARK_DATA_DIR) / name

    # Download file if not already cached or if force is True
    if force or not local_path.exists():
        print(f"[download] Downloading {name} from {uri}")
        response = requests.get(uri)
        response.raise_for_status()
        local_path.parent.mkdir(parents=True, exist_ok=True)
        with open(local_path, "wb") as f_out:
            f_out.write(response.content)
    else:
        print(f"[cache] Found cached file: {local_path.name}")

    # If it's a zip, extract it
    if local_path.suffix == ".zip":
        extract_dir = local_path.parent
        with zipfile.ZipFile(local_path, "r") as zip_ref:
            zip_ref.extractall(extract_dir)
        return extract_dir
    raise ValueError(f"Unsupported file format: {local_path}")


# OC157 benchmark node
class OC157Benchmark(zntrack.Node):
    """
    Benchmark model for OC157 dataset.

    Prediction of the most stable structures for a molecule-surface system
    - relative energies between 3 structures and 157 molecule surface combinations
    - identification of the most stable structure

    reference: MPRelaxSet DFT (Becke-Johnson damped D3 dispersion correction)
    - surfaces taken from the Open Catalyst Challenge 2023
    - 200 refs but excludes those with Hubbard U -> 157
    - 3 structures per system (triplet)
    """

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        """Run OC157 energy calculations."""
        calc = self.model.get_calculator()
        base_dir = get_benchmark_data("OC_Dataset.zip") / "OC_Dataset"
        n_systems = 200
        skip_hubbard_u = True

        # ---- Helper functions ----
        def _incar_has_hubbard_u(incar_path):
            if not incar_path.is_file():
                return False
            return "LDAU" in incar_path.read_text()

        def _find_energy(outcar, key="energy  without entropy="):
            with open(outcar, encoding="ISO-8859-1") as fh:
                hits = [line for line in fh if key in line]
            if not hits:
                raise RuntimeError(f"No energy found in {outcar}")
            return float(hits[-1].split()[-1])

        def _read_structure(folder):
            for fname in ("CONTCAR", "POSCAR"):
                fpath = folder / fname
                if fpath.is_file():
                    return read(fpath, format="vasp")
            raise FileNotFoundError(f"No CONTCAR/POSCAR in {folder}")

        def evaluate_energies(triplet):
            for atoms in triplet:
                atoms.calc = deepcopy(calc)
                atoms.get_potential_energy()

        triplets = []
        system_ids = []
        system_compositions = []

        for idx in tqdm(range(1, n_systems + 1), desc="Loading OC157 systems"):
            sys_id = f"{idx:03d}"
            sys_dir = base_dir / sys_id

            if skip_hubbard_u and _incar_has_hubbard_u(sys_dir / "1" / "INCAR"):
                continue

            poscar = (sys_dir / "1" / "POSCAR").read_text().splitlines()[0].strip()
            system_ids.append(sys_id)
            system_compositions.append(poscar)

            trio_atoms = []
            for member in (1, 2, 3):
                subdir = sys_dir / str(member)
                atoms = _read_structure(subdir)
                energy = _find_energy(subdir / "OUTCAR")
                atoms.info["ref_energy"] = energy
                atoms.info["composition"] = poscar
                atoms.info["sys_id"] = sys_id
                trio_atoms.append(atoms)
            triplets.append(trio_atoms)

        for trio in tqdm(triplets, desc="Evaluating model on triplets"):
            evaluate_energies(trio)
            write_dir = OUT_PATH / self.model_name
            write_dir.mkdir(parents=True, exist_ok=True)
            write(write_dir / f"{trio[-1].info['sys_id']}.xyz", trio)


project = mlipx.Project()
benchmark_node_dict = {}

for model_name, model in MODELS.items():
    with project.group(model_name):
        benchmark = OC157Benchmark(
            model=model,
            model_name=model_name,
        )
        benchmark_node_dict[model_name] = benchmark

project.build()
