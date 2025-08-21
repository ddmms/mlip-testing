"""Run S24 app."""

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

DATA_PATH = Path(__file__).parent.parent.parent / "data" / "surfaces" / "S24"

SCATTER = read_plot(DATA_PATH / "figure_adsorption_energies.json", id="figure")
STRUCTS_DIR = DATA_PATH / list(MODELS.keys())[0]
# Assets dir will be parent directory
STRUCTS = [
    f"assets/S24/{list(MODELS.keys())[0]}/{struct_file.stem}.xyz"
    for struct_file in STRUCTS_DIR.glob("*.xyz")
    if struct_file.name != "s24_mol_surface_atoms.xyz"
]


class S24App(BaseApp):
    """S24 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        plot_from_table_column(
            table_id=self.table_id,
            plot_id="figure-placeholder",
            column_to_plot={"MAE": SCATTER},
        )

        struct_from_scatter(
            scatter_id="figure",
            struct_id="struct-placeholder",
            structs=STRUCTS,
        )


def get_app() -> S24App:
    """
    Get S24 benchmark app layout and callback registration.

    Returns
    -------
    S24App
        Benchmark layout and callback registration.
    """
    return S24App(
        name="S24",
        title="S24",
        description=(
            "Performance in predicting adsorption energies for 24 "
            "molecule-surface combinations."
        ),
        table_path=DATA_PATH / "s24_metrics_table.json",
        extra_components=[Div(id="figure-placeholder"), Div(id="struct-placeholder")],
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
