"""Fixtures for MLIP results analysis."""

from __future__ import annotations

from collections.abc import Callable
import functools
from json import dump
from typing import Any

from dash import dash_table
import numpy as np
import plotly.graph_objects as go


def plot_parity(
    title: str | None = None,
    x_label: str | None = None,
    y_label: str | None = None,
    hoverdata: dict | None = (),
    filename: str = "parity.json",
) -> Callable:
    """
    Plot parity plot of MLIP results against reference data.

    Parameters
    ----------
    title
        Graph title.
    x_label
        Label for x-axis.
    y_label
        Label for y-axis.
    hoverdata
        Hover data dictionary.
    filename
        Filename to save table.

    Returns
    -------
    Callable
        Decorator to wrap function.
    """

    def plot_parity_decorator(func: Callable) -> Callable:
        """
        Decorate function to plot parity.

        Parameters
        ----------
        func
            Function being wrapped.

        Returns
        -------
        Callable
            Wrapped function.
        """

        @functools.wraps(func)
        def plot_parity_wrapper(*args, **kwargs) -> dict[str, Any]:
            """
            Wrap function to plot parity.

            Parameters
            ----------
            *args
                Arguments to pass to the function being wrapped.
            **kwargs
                Key word arguments to pass to the function being wrapped.

            Returns
            -------
            dict
                Results dictionary.
            """
            results = func(*args, **kwargs)
            ref = results["ref"]

            traces = []

            hovertemplate = "<b>Pred: </b>%{x}<br>" + "<b>Ref: </b>%{y}<br>"

            customdata = []
            for i, key in enumerate(hoverdata):
                hovertemplate += f"<b>{key}: </b>%{{customdata[{i}]}}<br>"
            customdata = list(zip(*hoverdata.values(), strict=True))

            for mlip, value in results.items():
                if mlip == "ref":
                    continue
                traces.append(
                    go.Scatter(
                        x=value,
                        y=ref,
                        name=mlip,
                        mode="markers",
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                    )
                )

            fig = go.Figure()

            for trace in traces:
                fig.add_trace(trace)

            full_fig = fig.full_figure_for_development()
            x_range = full_fig.layout.xaxis.range
            y_range = full_fig.layout.yaxis.range

            lims = [
                np.min([x_range, y_range]),  # min of both axes
                np.max([x_range, y_range]),  # max of both axes
            ]

            fig.add_trace(
                go.Scatter(
                    x=lims,
                    y=lims,
                    mode="lines",
                    showlegend=False,
                )
            )

            fig.update_layout(
                title={"text": title},
                xaxis={"title": {"text": x_label}},
                yaxis={"title": {"text": y_label}},
            )

            fig.update_traces()
            fig.write_json(filename)

            return results

        return plot_parity_wrapper

    return plot_parity_decorator


def build_table(
    filename: str = "table.json",
) -> Callable:
    """
    Build table MLIP results.

    Parameters
    ----------
    filename
        Filename to save table.

    Returns
    -------
    Callable
        Decorator to wrap function.
    """

    def build_table_decorator(func: Callable) -> Callable:
        """
        Decorate function to plot bar chart.

        Parameters
        ----------
        func
            Function being wrapped.

        Returns
        -------
        Callable
            Wrapped function.
        """

        @functools.wraps(func)
        def build_table_wrapper(*args, **kwargs) -> dict[str, Any]:
            """
            Wrap function to plot bar chart.

            Parameters
            ----------
            *args
                Arguments to pass to the function being wrapped.
            **kwargs
                Key word arguments to pass to the function being wrapped.

            Returns
            -------
            dict
                Results dictionary.
            """
            results = func(*args, **kwargs)
            # Form of results is
            # results = {
            #     metric_1: {mlip_1: value_1, mlip_2: value_2},
            #     metric_2: {mlip_1: value_3, mlip_2: value_4},
            # }

            metrics_columns = ("MLIP",) + tuple(results.keys())
            # Use MLIP keys from first (any) metric keys
            mlips = next(iter(results.values())).keys()

            metrics_data = []
            for mlip in mlips:
                metrics_data.append(
                    {"MLIP": mlip}
                    | {key: value[mlip] for key, value in results.items()},
                )

            table = dash_table.DataTable(
                metrics_data,
                [{"name": i, "id": i} for i in metrics_columns],
                id="metrics",
            )

            # Save dict of table to be loaded
            with open(filename, "w") as fp:
                dump({"data": table.data, "columns": table.columns}, fp)

            return results

        return build_table_wrapper

    return build_table_decorator
