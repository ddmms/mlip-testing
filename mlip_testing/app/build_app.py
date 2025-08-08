"""Build main Dash application."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path

from dash import Dash, Input, Output, callback
from dash.dash_table import DataTable
from dash.dcc import Store, Tab, Tabs
from dash.html import H1, Div

from mlip_testing import app
from mlip_testing.analysis.utils.utils import calc_ranks, calc_scores
from mlip_testing.app.utils.build_components import build_weight_components


def get_tabs() -> tuple[dict[str, list[Div]], dict[str, DataTable]]:
    """
    Get layout and register callbacks for all tab applications.

    Returns
    -------
    tuple[dict[str, list[Div]], dict[str, DataTable]]
        Layouts and tables for all tabs.
    """
    # Find Python files e.g. app_OC157.py in mlip_tesing.app module.
    tabs = Path(app.__file__).parent.glob("*/app*.py")
    layouts = {}
    tables = {}

    # Build all layouts, and register all callbacks to main app.
    for tab in tabs:
        # Import tab application layout/callbacks
        tab_name = tab.parent.name
        tab_module = import_module(f"mlip_testing.app.{tab_name}.app_{tab_name}")
        tab_app = tab_module.get_app()

        # Get layouts and tables for each tab
        layouts[tab_app.name] = tab_app.layout
        tables[tab_app.name] = tab_app.table

        # Register tab callbacks
        tab_app.register_callbacks()

    return layouts, tables


def build_summary_table(tables: dict[str, DataTable]) -> DataTable:
    """
    Build summary table.

    Parameters
    ----------
    tables
        Table from each tab.

    Returns
    -------
    DataTable
        Summary table with score from each tab.
    """
    summary_data = {}
    for tab_name, table in tables.items():
        summary_data = {row["MLIP"]: {} for row in table.data}
        for row in table.data:
            summary_data[row["MLIP"]][tab_name] = row["Score"]

    data = []
    for mlip in summary_data:
        data.append({"MLIP": mlip} | summary_data[mlip])

    data = calc_scores(data)
    data = calc_ranks(data)

    columns_headers = ("MLIP",) + tuple(tables.keys()) + ("Score", "Rank")
    columns = [{"name": headers, "id": headers} for headers in columns_headers]

    return DataTable(data=data, columns=columns, id="summary-table")


def build_tabs(
    full_app: Dash,
    layouts: dict[str, list[Div]],
    summary_table: DataTable,
    weight_components: Div,
) -> None:
    """
    Build tab layouts and summary tab.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.
    layouts
        Layouts for all tabs.
    summary_table
        Summary table with score from each tab.
    weight_components
        Weight sliders, text boxes and reset button.
    """
    all_tabs = [Tab(label="Summary", value="summary-tab", id="summary-tab")] + [
        Tab(label=tab_name, value=tab_name) for tab_name in layouts
    ]

    tabs_layout = [
        H1("MLIP benchmarking"),
        Tabs(id="all-tabs", value="summary-tab", children=all_tabs),
        Div(id="tabs-content"),
    ]

    full_app.layout = Div(tabs_layout)

    @callback(Output("tabs-content", "children"), Input("all-tabs", "value"))
    def select_tab(tab) -> Div:
        """
        Select tab contents to be displayed.

        Parameters
        ----------
        tab
            Name of tab selected.

        Returns
        -------
        Div
            Summary or tab contents to be displayed.
        """
        if tab == "summary-tab":
            return Div(
                [
                    H1("Benchmarks Summary"),
                    summary_table,
                    weight_components,
                    Store(
                        id="summary-table-scores-store",
                        storage_type="session",
                        # data=dict.fromkeys(columns, 1.0),
                    ),
                ]
            )
        return Div([layouts[tab]])


def build_full_app(full_app: Dash) -> None:
    """
    Build full app layout and register callbacks.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.
    """
    layouts, tables = get_tabs()
    summary_table = build_summary_table(tables)
    weight_components = build_weight_components(
        header="Benchmark weights",
        columns=list(tables.keys()),
        input_ids=list(tables.keys()),
        table_id="summary-table",
    )
    build_tabs(full_app, layouts, summary_table, weight_components)
