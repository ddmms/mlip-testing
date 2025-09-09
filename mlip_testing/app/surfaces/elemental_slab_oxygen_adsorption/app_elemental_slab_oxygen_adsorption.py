"""Run elemental_slab_oxygen_adsorption app."""

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

DATA_PATH = (
    Path(__file__).parent.parent.parent
    / "data"
    / "surfaces"
    / "elemental_slab_oxygen_adsorption"
)

SCATTER = read_plot(DATA_PATH / "figure_adsorption_energies.json", id="figure")
STRUCTS_DIR = DATA_PATH / list(MODELS.keys())[0]
# Assets dir will be parent directory
STRUCTS = [
    f"assets/elemental_slab_oxygen_adsorption/{list(MODELS.keys())[0]}/{struct_file.stem}.xyz"
    for struct_file in STRUCTS_DIR.glob("*.xyz")
]


class S24App(BaseApp):
    """S24 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        plot_from_table_column(
            table_id=self.table_id,
            plot_id="s24-figure-placeholder",
            column_to_plot={"MAE": SCATTER},
        )

        struct_from_scatter(
            scatter_id="figure",
            struct_id="s24-struct-placeholder",
            structs=STRUCTS,
        )


def get_app() -> S24App:
    """
    Get elemental_slab_oxygen_adsorption benchmark app layout and callback registration.

    Returns
    -------
    S24App
        Benchmark layout and callback registration.
    """
    return S24App(
        name="elemental_slab_oxygen_adsorption",
        title="Elemental Slab Oxygen Adsorption",
        description=(
            "Performance in predicting adsorption energies for oxygen on "
            "elemental slabs."
        ),
        table_path=DATA_PATH / "s24_metrics_table.json",
        extra_components=[
            Div(id="s24-figure-placeholder"),
            Div(id="s24-struct-placeholder"),
        ],
    )


if __name__ == "__main__":
    # Create Dash app
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent)

    # Construct layout and register callbacks
    s24_app = get_app()
    full_app.layout = s24_app.layout
    s24_app.register_callbacks()

    # Run app
    full_app.run(port=8052, debug=True)
