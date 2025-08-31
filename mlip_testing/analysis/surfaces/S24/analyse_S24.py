"""Analyse S24 benchmark."""

from __future__ import annotations

from pathlib import Path

from ase.io import read, write
import pytest

from mlip_testing.analysis.utils.decorators import build_table, plot_parity
from mlip_testing.analysis.utils.utils import mae
from mlip_testing.calcs.models.models import MODELS

CALC_PATH = (
    Path(__file__).parent.parent.parent.parent
    / "calcs"
    / "surfaces"
    / "S24"
    / "outputs"
)
OUT_PATH = (
    Path(__file__).parent.parent.parent.parent / "app" / "data" / "surfaces" / "S24"
)


def compute_adsorption_energy(
    surface_e: float, mol_surf_e: float, molecule_e: float
) -> float:
    """
    Compute adsorption energy.

    Parameters
    ----------
    surface_e
        Energy of the clean surface.
    mol_surf_e
        Energy of the molecule+surface system.
    molecule_e
        Energy of the isolated molecule.

    Returns
    -------
    float
        Adsorption energy.
    """
    return mol_surf_e - (surface_e + molecule_e)


def system_names() -> list:
    """
    Get list of system names.

    Returns
    -------
    list
        List of all system names.
    """
    system_names = []
    for model_name in MODELS:
        for system_path in (CALC_PATH / model_name).glob("*.xyz"):
            if system_path.name == "s24_mol_surface_atoms.xyz":
                continue
            system_names.append(system_path.stem)
        break
    return system_names


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_adsorption_energies.json",
    title="Adsorption energies",
    x_label="Predicted adsorption energy / eV",
    y_label="Reference adsorption energy / eV",
    hoverdata={
        "System": system_names(),
    },
)
def adsorption_energies() -> dict[str, list]:
    """
    Get adsorption energies for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted adsorption energies.
    """
    results = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False

    for model_name in MODELS:
        for system_path in sorted((CALC_PATH / model_name).glob("*.xyz")):
            if system_path.name == "s24_mol_surface_atoms.xyz":
                continue

            structs = read(system_path, index=":")
            if len(structs) != 3:
                continue

            surface, mol_surface, molecule = structs

            # Get predicted energies
            surface_e = surface.get_potential_energy()
            mol_surf_e = mol_surface.get_potential_energy()
            molecule_e = molecule.get_potential_energy()
            pred_ads_energy = compute_adsorption_energy(
                surface_e, mol_surf_e, molecule_e
            )
            results[model_name].append(pred_ads_energy)

            # Get reference energies (only store once)
            if not ref_stored:
                ref_surface_e = surface.info["ref_energy"]
                ref_mol_surf_e = mol_surface.info["ref_energy"]
                ref_molecule_e = molecule.info["ref_energy"]
                ref_ads_energy = compute_adsorption_energy(
                    ref_surface_e, ref_mol_surf_e, ref_molecule_e
                )
                results["ref"].append(ref_ads_energy)

            # Write structures in order
            structs_dir = OUT_PATH / model_name
            structs_dir.mkdir(parents=True, exist_ok=True)
            write(structs_dir / f"{system_path.stem}.xyz", structs)

        ref_stored = True
    return results


@pytest.fixture
def s24_mae(adsorption_energies) -> dict[str, float]:
    """
    Get mean absolute error for adsorption energies.

    Parameters
    ----------
    adsorption_energies
        Dictionary of reference and predicted adsorption energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted adsorption energy errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = mae(
            adsorption_energies["ref"], adsorption_energies[model_name]
        )
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "s24_metrics_table.json",
    metric_tooltips={
        "Model": "Name of the model",
        "MAE": "Mean Absolute Error (eV)",
    },
)
def metrics(s24_mae: dict[str, float]) -> dict[str, dict]:
    """
    Get all S24 metrics.

    Parameters
    ----------
    s24_mae
        Mean absolute errors for all models.

    Returns
    -------
    dict[str, dict]
        Metric names and values for all models.
    """
    return {
        "MAE": s24_mae,
    }


def test_s24(metrics: dict[str, dict]) -> None:
    """
    Run S24 test.

    Parameters
    ----------
    metrics
        All S24 metrics.
    """
    return
