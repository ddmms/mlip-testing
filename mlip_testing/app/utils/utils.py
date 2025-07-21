"""Utility functions for Dash app."""

from __future__ import annotations

import json
from pathlib import Path

from dash.dash_table import DataTable
from plotly.graph_objs import Figure
from plotly.io import read_json


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


def read_plot(filename: str | Path) -> Figure:
    """
    Read preprepared plotly Figure.

    Parameters
    ----------
    filename
        Name of json file with saved plot data.

    Returns
    -------
    Figure
        Loaded plotly Figure.
    """
    return read_json(filename)
