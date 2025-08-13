"""Helpers to create callbaclks for Dash app."""

from __future__ import annotations

from dash import Input, Output, callback
from dash.dcc import Graph
from dash.html import Div, Iframe

from mlip_testing.app.utils.weas import generate_weas_html


def plot_from_table_column(
    table_id: str, plot_id: str, column_to_plot: dict[str, Graph]
) -> None:
    """
    Attach callback to show plot when a table column is clicked.

    Parameters
    ----------
    table_id
        ID for Dash table being clicked.
    plot_id
        ID for Dash plot placeholder Div.
    column_to_plot
        Dictionary relating table headers (keys) and plot to show (values).
    """

    @callback(Output(plot_id, "children"), Input(table_id, "active_cell"))
    def show_plot(active_cell) -> Div:
        """
        Register callback to show plot when a table column is clicked.

        Parameters
        ----------
        active_cell
            Clicked cell in Dash table.

        Returns
        -------
        Div
            Message explaining interactivity, or plot on table click.
        """
        if not active_cell:
            return Div("Click on a metric to view plot.")
        column_id = active_cell.get("column_id", None)
        if column_id:
            return Div(column_to_plot[column_id])
        raise ValueError("Invalid column_id")


def plot_from_table_cell(
    table_id: str,
    plot_id: str,
    cell_to_plot: dict[str, dict[Graph]],
) -> None:
    """
    Attach callback to show plot when a table cell is clicked.

    Parameters
    ----------
    table_id
        ID for Dash table being clicked.
    plot_id
        ID for Dash plot placeholder Div.
    cell_to_plot
        Nested dictionary of model names, column names, and plot to show.
    """

    @callback(Output(plot_id, "children"), Input(table_id, "active_cell"))
    def show_plot(active_cell) -> Div:
        """
        Register callback to show plot when a table cell is clicked.

        Parameters
        ----------
        active_cell
            Clicked cell in Dash table.

        Returns
        -------
        Div
            Message explaining interactivity, or plot on cell click.
        """
        if not active_cell:
            return Div("Click on a metric to view plot.")
        column_id = active_cell.get("column_id", None)
        row_id = active_cell.get("row_id", None)

        if row_id in cell_to_plot and column_id in cell_to_plot[row_id]:
            return Div(cell_to_plot[row_id][column_id])
        return Div("Click on a metric to view plot.")


def struct_from_scatter(scatter_id: str, struct_id: str, structs: list[str]) -> None:
    """
    Attach callback to show a structure when a scatter point is clicked.

    Parameters
    ----------
    scatter_id
        ID for Dash scatter being clicked.
    struct_id
        ID for Dash plot placeholder Div.
    structs
        List of structure filenames in same order as scatter data to be visualised.
    """

    @callback(Output(struct_id, "children"), Input(scatter_id, "clickData"))
    def show_struct(click_data):
        """
        Register callback to show structure when a scatter point is clicked.

        Parameters
        ----------
        click_data
            Clicked data point in scatter plot.

        Returns
        -------
        Div
            Visualised structure on plot click.
        """
        if not click_data:
            return None
        idx = click_data["points"][0]["pointNumber"]
        return Div(
            Iframe(
                srcDoc=generate_weas_html(structs[idx]),
                style={
                    "height": "550px",
                    "width": "100%",
                    "border": "1px solid #ddd",
                    "borderRadius": "5px",
                },
            )
        )
