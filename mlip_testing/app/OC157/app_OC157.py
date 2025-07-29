"""Run OC157 app."""

from __future__ import annotations

from pathlib import Path

from dash import Dash
from dash.html import Div

from mlip_testing.app.utils.build import layout_builder
from mlip_testing.app.utils.callbacks import plot_from_table_column, struct_from_scatter
from mlip_testing.app.utils.load import read_plot, rebuild_table
from mlip_testing.calcs.models.models import MODELS

DATA_PATH = Path(__file__).parent.parent / "data" / "OC157"


metrics_table = rebuild_table(DATA_PATH / "oc157_metrics_table.json", id="table")
scatter = read_plot(DATA_PATH / "figure_rel_energies.json", id="figure")
structs_dir = DATA_PATH / list(MODELS.keys())[0]
# Assets dir will be parent directory
structs = [
    f"assets/OC157/{list(MODELS.keys())[0]}/{i}.xyz"
    for i in range(len(list(structs_dir.glob("*.xyz"))))
]


def build_layout() -> Div:
    """
    Build app layout.

    Returns
    -------
    Div
        Div component with list all components for app.
    """
    # Define all components/placeholders
    return layout_builder(
        title="OC157 benchmark",
        table=metrics_table,
        extra_components=[Div(id="figure-placeholder"), Div(id="struct-placeholder")],
    )


def register_callbacks(app: Dash) -> None:
    """
    Register callbacks to app.

    Parameters
    ----------
    app
        Dash application to register callbacks to.
    """
    plot_from_table_column(
        app=app,
        table_id="table",
        plot_id="figure-placeholder",
        column_to_plot={"MAE (meV)": scatter, "Ranking Error": scatter},
    )

    struct_from_scatter(
        app=app,
        scatter_id="figure",
        struct_id="struct-placeholder",
        structs=structs,
    )


if __name__ == "__main__":
    # Create and run Dash app
    app = Dash(__name__, assets_folder=DATA_PATH.parent)
    app.layout = build_layout()
    register_callbacks(app)
    app.run(port=8051, debug=True)
