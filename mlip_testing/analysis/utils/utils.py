"""Utility functions for analysis."""

from __future__ import annotations

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
