"""Run OC157 app."""

from __future__ import annotations

from pathlib import Path

from dash import Dash
from dash.html import Div

from mlip_testing.app.base_app import BaseApp
from mlip_testing.app.utils.build_callbacks import (
    plot_from_table_column,
    struct_from_scatter,
)
from mlip_testing.app.utils.load import read_plot
from mlip_testing.calcs.models.models import MODELS

DATA_PATH = Path(__file__).parent.parent / "data" / "OC157"

SCATTER = read_plot(DATA_PATH / "figure_rel_energies.json", id="figure")
STRUCTS_DIR = DATA_PATH / list(MODELS.keys())[0]
# Assets dir will be parent directory
STRUCTS = [
    f"assets/OC157/{list(MODELS.keys())[0]}/{i}.xyz"
    for i in range(len(list(STRUCTS_DIR.glob("*.xyz"))))
]


class OC157App(BaseApp):
    """OC157 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        plot_from_table_column(
            table_id=self.table_id,
            plot_id="figure-placeholder",
            column_to_plot={"MAE": SCATTER, "Ranking Error": SCATTER},
        )

        struct_from_scatter(
            scatter_id="figure",
            struct_id="struct-placeholder",
            structs=STRUCTS,
        )


def get_app() -> OC157App:
    """
    Get OC157 benchmark app layout and callback registration.

    Returns
    -------
    OC157App
        Benchmark layout and callback registration.
    """
    return OC157App(
        name="OC157",
        title="OC157",
        description=(
            "Performance in predicting relative energies between 3 structures for 157 "
            "molecule-surface combinations."
        ),
        table_path=DATA_PATH / "oc157_metrics_table.json",
        extra_components=[Div(id="figure-placeholder"), Div(id="struct-placeholder")],
    )


if __name__ == "__main__":
    # Create Dash app
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent)

    # Construct layout and register callbacks
    oc157_app = get_app()
    full_app.layout = oc157_app.layout
    oc157_app.register_callbacks()

    # Run app
    full_app.run(port=8051, debug=True)
