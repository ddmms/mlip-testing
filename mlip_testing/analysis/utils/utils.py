"""Utility functions for analysis."""

from __future__ import annotations

import json
from pathlib import Path

from dash.dash_table import DataTable
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


def rebuild_table(filename: str | Path) -> DataTable:
    """
    Rebuild saved dash table.

    Parameters
    ----------
    filename
        Name of json file with saved table data.

    Returns
    -------
    DataTable
        Dash DataTable.
    """
    # Load JSON file
    with open(filename) as f:
        table_json = json.load(f)

    data = table_json["data"]
    columns = table_json["columns"]

    return DataTable(data=data, columns=columns, editable=True)
