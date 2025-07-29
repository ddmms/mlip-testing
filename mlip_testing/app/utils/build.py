"""Utility functions for building app components."""

from __future__ import annotations

from dash.dash_table import DataTable
from dash.development.base_component import Component
from dash.html import H1, Div


def layout_builder(
    title: str, table: DataTable, extra_components: list[Component] | None = None
) -> Div:
    """
    Build app layout.

    Parameters
    ----------
    title
        Title for app tab.
    table
        Dash Table with metric results.
    extra_components
        List of Dash Components to include after the metrics table.

    Returns
    -------
    Div
        App tab layout.
    """
    layout_contents = [H1(title, style={"color": "black"}), Div(table)]
    if extra_components:
        layout_contents.extend(extra_components)

    return Div(layout_contents)
