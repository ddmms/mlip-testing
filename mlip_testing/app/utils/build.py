"""Utility functions for building app components."""

from __future__ import annotations

from dash import Input, Output, callback, ctx
from dash.dash_table import DataTable
from dash.dcc import Input as DCC_Input
from dash.dcc import Slider, Store
from dash.development.base_component import Component
from dash.exceptions import PreventUpdate
from dash.html import H1, H2, Button, Div, Label

from mlip_testing.analysis.utils.utils import calc_ranks, calc_scores


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


def register_table_callback(table_id) -> None:
    """
    Register callback to update table data.

    Parameters
    ----------
    table_id
        ID for table to update.
    """

    # Callback update table when store changes.
    @callback(
        Output(f"{table_id}", "data"),
        Input(f"{table_id}-weight-store", "data"),
        Input(f"{table_id}", "data"),
        prevent_initial_call=False,
    )
    def update_table(stored_values: dict[str:float], table_data: list[dict]) -> dict:
        """
        Update table contents from data store.

        Parameters
        ----------
        stored_values
            Stored value.
        table_data
            Data from table to be updated.

        Returns
        -------
        list[dict]
            Updated table data.
        """
        trigger_id = ctx.triggered_id

        if trigger_id != f"{table_id}-weight-store" and trigger_id is not None:
            raise ValueError("Invalid trigger. trigger_id: ", trigger_id)

        # Update table contents
        table_data = calc_scores(table_data, stored_values)
        return calc_ranks(table_data)


def register_weight_callbacks(input_id: str, table_id: str, column: str) -> None:
    """
    Register all callbacks for weight inputs.

    Parameters
    ----------
    input_id
        ID prefix for slider and input box.
    table_id
        ID for table. Also used to identify reset button and weight store.
    column
        Column header corresponding to slider and input box.
    """
    default_value = 1.0

    @callback(
        Output(
            f"{table_id}-weight-store",
            "data",
            allow_duplicate=True,
        ),
        Input(f"{input_id}-slider", "value"),
        Input(f"{table_id}-weight-store", "data"),
        prevent_initial_call=True,
    )
    def store_slider_value(
        slider_value: float, stored_values: dict[str:float]
    ) -> dict[str, float]:
        """
        Store weight values from slider and text input.

        Parameters
        ----------
        slider_value
            Value from slider.
        stored_values
            Stored values dictionary.

        Returns
        -------
        dict[str, float]
            Stored weights for each slider.
        """
        trigger_id = ctx.triggered_id

        if trigger_id == f"{input_id}-slider":
            weight = slider_value
        else:
            raise PreventUpdate
        stored_values[column] = weight

        return stored_values

    @callback(
        Output(f"{input_id}-slider", "value", allow_duplicate=True),
        Output(f"{input_id}-input", "value", allow_duplicate=True),
        Input(f"{table_id}-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_slider(n_clicks) -> tuple[float, float]:
        """
        Reset slider value on click.

        Parameters
        ----------
        n_clicks
            Number of clicks of reset button.

        Returns
        -------
        tuple[float, float]
            Weight to set slider and input box value.
        """
        return default_value, default_value

    @callback(
        Output(f"{input_id}-slider", "value", allow_duplicate=True),
        Output(f"{input_id}-input", "value", allow_duplicate=True),
        Input(f"{table_id}-weight-store", "data"),
        Input("all-tabs", "value"),
        prevent_initial_call="initial_duplicate",
    )
    def load_store(
        stored_values: dict[str, float], tabs_value: str
    ) -> tuple[float, float]:
        """
        Load stored values.

        Parameters
        ----------
        stored_values
            Number of clicks of reset button.
        tabs_value
            Tab name.

        Returns
        -------
        tuple[float, float]
            Weight to set slider and input box value.
        """
        return stored_values[column], stored_values[column]

    @callback(
        Output(f"{input_id}-input", "value", allow_duplicate=True),
        Output(f"{input_id}-slider", "value", allow_duplicate=True),
        Input(f"{input_id}-input", "value"),
        Input(f"{input_id}-slider", "value"),
        prevent_initial_call=True,
    )
    def sync_slider_input(input_value, slider_value) -> tuple[float, float, float]:
        """
        Sync weight values between slider and text input.

        Parameters
        ----------
        input_value
            Value from text box input.
        slider_value
            Value from slider.

        Returns
        -------
        tuple[float, float]
            Weights to set slider value and text input value.
        """
        trigger_id = ctx.triggered_id

        if trigger_id == f"{input_id}-input":
            if input_value:
                weight = input_value
            else:
                raise PreventUpdate
        elif trigger_id == f"{input_id}-slider":
            weight = slider_value
        else:
            raise ValueError("Invalid trigger. trigger_id: ", trigger_id)

        return weight, weight


def build_weight_components(
    header: str,
    columns: list[str],
    input_ids: list[str],
    table_id: str,
) -> Div:
    """
    Build weight sliders, text boxes and reset button.

    Parameters
    ----------
    header
        Header for above sliders.
    columns
        Column headers to look up stored values, and label sliders and input boxes.
    input_ids
        ID prefixes for sliders and input boxes.
    table_id
        ID for table. Also used to identify reset button and weight store.

    Returns
    -------
    Div
        Div containing header, weight sliders, text boxes and reset button.
    """
    layout = [Div(header)]

    for column, input_id in zip(columns, input_ids, strict=True):
        layout.append(
            build_slider(
                label=column,
                slider_id=f"{input_id}-slider",
                input_id=f"{input_id}-input",
                default_value=None,  # Set by stored value/default
            )
        )

    layout.extend(
        [
            Button(
                "Reset Weights",
                id=f"{table_id}-reset-button",
                n_clicks=0,
                style={"marginTop": "20px"},
            ),
            Store(
                id=f"{table_id}-weight-store",
                storage_type="session",
                data=dict.fromkeys(columns, 1.0),
            ),
        ]
    )

    register_table_callback(table_id=table_id)
    for column, input_id in zip(columns, input_ids, strict=True):
        register_weight_callbacks(input_id=input_id, table_id=table_id, column=column)

    return Div(layout)


def build_tab(
    name: str,
    title: str,
    description: str,
    table: DataTable,
    table_id: str,
    extra_components: list[Component] | None = None,
) -> Div:
    """
    Build app tab layout.

    Parameters
    ----------
    name
        Name for application tab.
    title
        Title for app tab.
    description
        Description of benchmark.
    table
        Dash Table with metric results.
    table_id
        ID of Dash Table.
    extra_components
        List of Dash Components to include after the metrics table.

    Returns
    -------
    Div
        App tab layout.
    """
    layout_contents = [
        H1(title, style={"color": "black"}),
        H2(description),
        Div(table),
    ]

    metrics = [
        col["name"]
        for col in table.columns
        if col["name"] not in ("MLIP", "Score", "Rank")
    ]
    ids = [f"{name}-{metric}" for metric in metrics]

    layout_contents.append(
        build_weight_components(
            header="Metric Weights",
            columns=metrics,
            input_ids=ids,
            table_id=table_id,
        )
    )

    if extra_components:
        layout_contents.extend(extra_components)

    return Div(layout_contents)
