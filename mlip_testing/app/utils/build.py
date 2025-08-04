"""Utility functions for building app components."""

from __future__ import annotations

from dash import Input, Output, callback, ctx
from dash.dash_table import DataTable
from dash.dcc import Input as DCC_Input
from dash.dcc import Slider, Store
from dash.development.base_component import Component
from dash.html import H1, H2, Button, Div, Label


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


def register_weight_callbacks(weight_prefix: str, table_prefix: str) -> None:
    """
    Register all callbacks for weight inputs.

    Parameters
    ----------
    weight_prefix
        Prefix for weight component IDs.
    table_prefix
        Prefix for table and reset component ID.
    """

    # Callback to sync weights between slider, text, reset, and Store
    @callback(
        Output(f"{weight_prefix}-input", "value"),
        Output(f"{weight_prefix}-slider", "value"),
        Output(f"{weight_prefix}-weight-store", "data"),
        Input(f"{weight_prefix}-input", "value"),
        Input(f"{weight_prefix}-slider", "value"),
        Input(f"{weight_prefix}-weight-store", "data"),
        Input(f"{table_prefix}-reset-weights-button", "n_clicks"),
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

        if trigger_id == f"{weight_prefix}-weight-store" or trigger_id is None:
            weight = store
        elif trigger_id == f"{weight_prefix}-input":
            weight = input_value
        elif trigger_id == f"{weight_prefix}-slider":
            weight = slider_value
        elif trigger_id == f"{table_prefix}-reset-weights-button":
            weight = 1
        else:
            raise ValueError("Invalid trigger. trigger_id: ", trigger_id)
        return weight, weight, weight


def build_weight_components(
    header: str, labels: list[str], ids: list[str], table_prefix: str
) -> Div:
    """
    Build weight sliders, text boxes and reset button.

    Parameters
    ----------
    header
        Header for above sliders.
    labels
        Names for each weight slider.
    ids
        Prefix for slider and input box IDs.
    table_prefix
        Label for table and reset button used for all weights.

    Returns
    -------
    Div
        Div containing header, weight sliders, text boxes and reset button.
    """
    layout = [Div(header)]

    for label, weight_prefix in zip(labels, ids, strict=True):
        layout.append(
            build_slider(
                label=label,
                slider_id=f"{weight_prefix}-slider",
                input_id=f"{weight_prefix}-input",
                default_value=None,  # Set by stored value/default
            )
        )

    layout.append(
        Button(
            "Reset Weights",
            id=f"{table_prefix}-reset-weights-button",
            n_clicks=0,
            style={"marginTop": "20px"},
        ),
    )

    for weight_prefix in ids:
        layout.append(
            Store(id=f"{weight_prefix}-weight-store", storage_type="session", data=1.0)
        )
        register_weight_callbacks(weight_prefix, table_prefix)

    return Div(layout)


def build_tab(
    name: str,
    title: str,
    description: str,
    table: DataTable,
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
            labels=metrics,
            ids=ids,
            reset_prefix=name,
        )
    )

    if extra_components:
        layout_contents.extend(extra_components)

    return Div(layout_contents)
