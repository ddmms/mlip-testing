"""Analyse Li diffusion benchmark."""

from __future__ import annotations

from pathlib import Path

from ase.io import read, write
import pytest

from mlip_testing.analysis.utils.decorators import build_table
from mlip_testing.calcs.models.models import MODELS

CALC_PATH = (
    Path(__file__).parent.parent.parent.parent
    / "calcs"
    / "nebs"
    / "li_diffusion"
    / "outputs"
)
OUT_PATH = Path(__file__).parent.parent.parent.parent / "app" / "data" / "li_diffusion"

REF_VALUES = {"path_b": 0.27, "path_c": 2.5}


@pytest.fixture
def path_b_error() -> dict[str, float]:
    """
    Get error in path B energy barrier.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted barrier errors for all models.
    """
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    results = {}
    for model_name in MODELS:
        structs = read(
            CALC_PATH / f"li_diffusion_b-{model_name}-neb-band.extxyz", index=":"
        )
        write(OUT_PATH / f"{model_name}-b-neb-band.extxyz", structs)

        with open(
            CALC_PATH / f"li_diffusion_b-{model_name}-neb-results.dat", encoding="utf8"
        ) as f:
            data = f.readlines()
            pred_barrier, _, _ = tuple(float(x) for x in data[1].split())
        results[model_name] = abs(REF_VALUES["path_b"] - pred_barrier)
    return results


@pytest.fixture
def path_c_error() -> dict[str, float]:
    """
    Get error in path B energy barrier.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted barrier errors for all models.
    """
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    results = {}
    for model_name in MODELS:
        structs = read(
            CALC_PATH / f"li_diffusion_c-{model_name}-neb-band.extxyz", index=":"
        )
        write(OUT_PATH / f"{model_name}-c-neb-band.extxyz", structs)

        with open(
            CALC_PATH / f"li_diffusion_c-{model_name}-neb-results.dat", encoding="utf8"
        ) as f:
            data = f.readlines()
            pred_barrier, _, _ = tuple(float(x) for x in data[1].split())
        results[model_name] = abs(REF_VALUES["path_c"] - pred_barrier)
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "li_diffusion_metrics_table.json",
    metric_tooltips={
        "Model": "Name of the model",
        "Path B error": "Energy Barrier error for path B (eV)",
        "Path C error": "Energy Barrier error for path C (eV)",
    },
)
def metrics(
    path_b_error: dict[str, float], path_c_error: dict[str, float]
) -> dict[str, dict]:
    """
    Get all Li diffusion metrics.

    Parameters
    ----------
    path_b_error
        Mean absolute errors for all models.
    path_c_error
        Mean absolute errors for all models.

    Returns
    -------
    dict[str, dict]
        Metric names and values for all models.
    """
    return {
        "Path B error": path_b_error,
        "Path C error": path_c_error,
    }


def test_li_diffusion(metrics: dict[str, dict]) -> None:
    """
    Run Li diffusion test.

    Parameters
    ----------
    metrics
        All Li diffusion metrics.
    """
    return
