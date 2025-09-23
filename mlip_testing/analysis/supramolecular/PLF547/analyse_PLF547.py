"""Analyse PLF547 benchmark."""

from __future__ import annotations

from ase import units
from ase.io import read, write
import numpy as np
import pytest

from mlip_testing.analysis.utils.decorators import build_table, plot_parity
from mlip_testing.analysis.utils.utils import mae
from mlip_testing.app import APP_ROOT
from mlip_testing.calcs import CALCS_ROOT
from mlip_testing.calcs.models.models import MODELS

CALC_PATH = CALCS_ROOT / "supramolecular" / "PLF547" / "outputs"
OUT_PATH = APP_ROOT / "data" / "supramolecular" / "PLF547"

# Constants
KCAL_PER_MOL_TO_EV = units.kcal / units.mol
EV_TO_KCAL_PER_MOL = 1.0 / KCAL_PER_MOL_TO_EV


def get_interaction_energies() -> dict[str, list]:
    """
    Get interaction energies for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted interaction energies.
    """
    results = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False

    for model_name in MODELS:
        model_dir = CALC_PATH / model_name
        if not model_dir.exists():
            continue

        # Sort files to ensure consistent ordering across models
        system_paths = sorted(model_dir.glob("*.xyz"))
        for system_path in system_paths:
            try:
                atoms = read(system_path)
                pred_energy = atoms.info["interaction_energy"]
                results[model_name].append(pred_energy)

                if not ref_stored:
                    ref_energy = atoms.info["ref_interaction_energy"]
                    results["ref"].append(ref_energy)

            except (KeyError, FileNotFoundError):
                continue

        ref_stored = True

    return results


def compositions() -> list:
    """
    Get list of system identifiers.

    Returns
    -------
    list
        List of all system identifiers.
    """
    all_compositions = []
    for model_name in MODELS:
        model_dir = CALC_PATH / model_name
        if not model_dir.exists():
            continue

        # Sort files to ensure consistent ordering
        system_paths = sorted(model_dir.glob("*.xyz"))
        for system_path in system_paths:
            try:
                atoms = read(system_path)
                sys_id = atoms.info.get("sys_id", system_path.stem)
                all_compositions.append(sys_id)
            except (OSError, KeyError, ValueError):
                all_compositions.append(system_path.stem)
        break
    return all_compositions


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_interaction_energies.json",
    title="Interaction energies",
    x_label="Predicted interaction energy / kcal/mol",
    y_label="Reference interaction energy / kcal/mol",
    hoverdata={
        "System": compositions(),
    },
)
def interaction_energies() -> dict[str, list]:
    """
    Get interaction energies for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted interaction energies in kcal/mol.
    """
    energies = get_interaction_energies()

    # Convert to kcal/mol
    results = {}
    for key, values in energies.items():
        results[key] = [e * EV_TO_KCAL_PER_MOL for e in values]

    # Copy structures to output directory for visualization
    for model_name in MODELS:
        model_dir = CALC_PATH / model_name
        if not model_dir.exists():
            continue

        structs_dir = OUT_PATH / model_name
        structs_dir.mkdir(parents=True, exist_ok=True)

        # Sort files to ensure consistent ordering
        system_paths = sorted(model_dir.glob("*.xyz"))
        for i, system_path in enumerate(system_paths):
            atoms = read(system_path)
            write(structs_dir / f"{i}.xyz", atoms)
        break  # Only need to copy once

    return results


@pytest.fixture
def plf547_mae(interaction_energies) -> dict[str, float]:
    """
    Get mean absolute error across all systems.

    Parameters
    ----------
    interaction_energies
        Dictionary of reference and predicted interaction energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted interaction energy MAE for all models in kcal/mol.
    """
    results = {}
    for model_name in MODELS:
        if (
            model_name in interaction_energies
            and len(interaction_energies[model_name]) > 0
        ):
            results[model_name] = mae(
                interaction_energies["ref"], interaction_energies[model_name]
            )
        else:
            results[model_name] = np.nan
    return results


@pytest.fixture
def r_squared(interaction_energies) -> dict[str, float]:
    """
    Get R-squared (coefficient of determination) across all systems.

    Parameters
    ----------
    interaction_energies
        Dictionary of reference and predicted interaction energies.

    Returns
    -------
    dict[str, float]
        Dictionary of R² values for all models.
    """
    results = {}
    for model_name in MODELS:
        if (
            model_name in interaction_energies
            and len(interaction_energies[model_name]) > 0
        ):
            ref = np.array(interaction_energies["ref"])
            pred = np.array(interaction_energies[model_name])

            # Calculate R²
            ss_res = np.sum((ref - pred) ** 2)
            ss_tot = np.sum((ref - np.mean(ref)) ** 2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0

            results[model_name] = r2
        else:
            results[model_name] = 0.0
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "plf547_metrics_table.json",
    metric_tooltips={
        "Model": "Name of the model",
        "MAE": "Mean Absolute Error (kcal/mol)",
        "R²": "Coefficient of determination",
    },
)
def metrics(
    plf547_mae: dict[str, float],
    r_squared: dict[str, float],
) -> dict[str, dict]:
    """
    Get all PLF547 metrics.

    Parameters
    ----------
    plf547_mae
        Mean absolute errors for all models.
    r_squared
        R-squared values for all models.

    Returns
    -------
    dict[str, dict]
        Metric names and values for all models.
    """
    return {
        "MAE": plf547_mae,
        "R²": r_squared,
    }


def test_plf547(metrics: dict[str, dict]) -> None:
    """
    Run PLF547 test.

    Parameters
    ----------
    metrics
        All PLF547 metrics.
    """
    return
