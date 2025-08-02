"""Build main Dash application."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path

from dash import Dash, Input, Output, callback, ctx
from dash.dash_table import DataTable
from dash.dcc import Input as DCC_Input
from dash.dcc import Slider, Store, Tab, Tabs
from dash.html import H1, Button, Div, Label

from mlip_testing import app
from mlip_testing.analysis.utils.utils import calc_ranks, calc_scores


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
        layouts[tab_name] = tab_app.layout
        tables[tab_name] = tab_app.table

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

    return DataTable(data=data, columns=columns)


def build_slider(
    label: str, slider_id: str, input_id: str, default_value: float | None
) -> Div:
    """
    Build slider and input box.

    Parameters
    ----------
    label
        Slider label.
    slider_id
        ID for slider component.
    input_id
        ID for text box input component.
    default_value
        Default value for slider/text box input.

    Returns
    -------
    Div
        Slider and input text box.
    """
    return Div(
        [
            Label(label),
            Div(
                [
                    Div(
                        Slider(
                            id=slider_id,
                            min=0,
                            max=5,
                            step=0.1,
                            value=default_value,
                            tooltip={"always_visible": False},
                            marks=None,
                        ),
                        style={"flex": "1 1 80%"},
                    ),
                    DCC_Input(
                        id=input_id,
                        type="number",
                        value=default_value,
                        step=0.1,
                        style={"width": "80px"},
                    ),
                ],
                style={"display": "flex", "gap": "10px", "alignItems": "center"},
            ),
        ]
    )


def register_weight_callbacks(tab_name: str) -> None:
    """
    Register all callbacks for weight inputs.

    Parameters
    ----------
    tab_name
        Name of tab.
    """

    # Callback to sync weights between slider, text, reset, and Store
    @callback(
        Output(f"{tab_name}-input", "value"),
        Output(f"{tab_name}-slider", "value"),
        Output(f"{tab_name}-weight-store", "data"),
        Input(f"{tab_name}-input", "value"),
        Input(f"{tab_name}-slider", "value"),
        Input(f"{tab_name}-weight-store", "data"),
        Input("reset-weights-button", "n_clicks"),
        prevent_initial_call=False,
    )
    def store_slider_value(
        input_value, slider_value, store, reset
    ) -> tuple[float, float, float]:
        """
        Store, reset, and sync weight values between slider and text input.

        Parameters
        ----------
        input_value
            Value from text box input.
        slider_value
            Value from slider.
        store
            Stored value.
        reset
            Number of clicks of reset button.

        Returns
        -------
        tuple[float, float, float]
            Weights to set slider value, text input value and stored value.
        """
        trigger_id = ctx.triggered_id

        if trigger_id == f"{tab_name}-weight-store":
            weight = store
        elif trigger_id == f"{tab_name}-input":
            weight = input_value
        elif trigger_id == f"{tab_name}-slider":
            weight = slider_value
        elif trigger_id == "reset-weights-button":
            weight = 1
        else:
            raise ValueError("Invalid trigger. trigger_id: ", trigger_id)
        return weight, weight, weight


def build_weight_components(tables: dict[str, DataTable]) -> Div:
    """
    Build weight sliders, text boxes and reset button.

    Parameters
    ----------
    tables
        Table from each tab.

    Returns
    -------
    Div
        Div containing weight sliders, text boxes and reset button.
    """
    layout = [Div("Benchmark Weights")]

    for tab_name in tables:
        layout.append(
            build_slider(
                label=tab_name,
                slider_id=f"{tab_name}-slider",
                input_id=f"{tab_name}-input",
                default_value=None,  # Set by stored value/default
            )
        )

    layout.append(
        Button(
            "Reset Weights",
            id="reset-weights-button",
            n_clicks=0,
            style={"marginTop": "20px"},
        ),
    )

    layout.append(
        Store(id=f"{tab_name}-weight-store", storage_type="session", data=1.0)
    )

    for tab_name in tables:
        register_weight_callbacks(tab_name)

    return Div(layout)


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
    all_tabs = [Tab(label="Summary", value="summary-tab")] + [
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
    weight_components = build_weight_components(tables)
    build_tabs(full_app, layouts, summary_table, weight_components)
