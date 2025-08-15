"""Utility functions for analysis."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error


def mae(ref: list, prediction: list) -> float:
    """
    Get mean absolute error.

    Parameters
    ----------
    ref
        Reference data.
    prediction
        Predicted data.

    Returns
    -------
    float
        Mean absolute error.
    """
    return mean_absolute_error(ref, prediction)


def rmse(ref: list, prediction: list) -> float:
    """
    Get root mean squared error.

    Parameters
    ----------
    ref
        Reference data.
    prediction
        Predicted data.

    Returns
    -------
    float
        Root mean squared error.
    """
    return mean_squared_error(ref, prediction)


def calc_scores(
    metrics_data: list[dict], weights: dict[str, float] | None = None
) -> list[dict]:
    """
    Calculate score for each model and add to table data.

    Parameters
    ----------
    metrics_data
        Rows data containing model name and metric values.
    weights
        Weight for each metric. Default is 1.0 for each metric.

    Returns
    -------
    list[dict]
        Rows of data with combined score for each model added.
    """
    weights = weights if weights else {}

    for row in metrics_data:
        score = 0
        for key, value in row.items():
            weight = weights.get(key, 1.0)
            if key in ("MLIP", "Score", "Rank"):
                continue
            score += value * weight
        row["Score"] = score

    return metrics_data


def calc_ranks(metrics_data: list[dict]) -> list[dict]:
    """
    Calculate rank for each model and add to table data.

    Parameters
    ----------
    metrics_data
        Rows data containing model name, metric values, and Score.

    Returns
    -------
    list[dict]
        Rows of data with rank for each model added.
    """
    ranked_scores = np.argsort([x["Score"] for x in metrics_data]) + 1
    for i, row in enumerate(metrics_data):
        row["Rank"] = ranked_scores[i]
    return metrics_data
