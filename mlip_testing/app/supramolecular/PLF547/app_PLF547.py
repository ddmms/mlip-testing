"""Run PLF547 app."""

from __future__ import annotations

from pathlib import Path

from dash import Dash
from dash.html import Div

from mlip_testing.app import APP_ROOT
from mlip_testing.app.base_app import BaseApp
from mlip_testing.app.utils.build_callbacks import (
    plot_from_table_column,
    struct_from_scatter,
)
from mlip_testing.app.utils.load import read_plot
from mlip_testing.calcs.models.models import MODELS

BENCHMARK_NAME = Path(__file__).name.removeprefix("app_").removesuffix(".py")
DATA_PATH = APP_ROOT / "data" / "supramolecular" / "PLF547"


class PLF547App(BaseApp):
    """PLF547 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        scatter = read_plot(
            DATA_PATH / "figure_interaction_energies.json",
            id=f"{BENCHMARK_NAME}-figure",
        )

        structs_dir = DATA_PATH / list(MODELS.keys())[0]

        # Get structure paths for visualization
        structs = []
        if structs_dir.exists():
            structs = [
                f"assets/supramolecular/PLF547/{list(MODELS.keys())[0]}/{i}.xyz"
                for i in range(len(list(structs_dir.glob("*.xyz"))))
            ]

        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={
                "MAE": scatter,
                "R²": scatter,
            },
        )

        struct_from_scatter(
            scatter_id=f"{BENCHMARK_NAME}-figure",
            struct_id=f"{BENCHMARK_NAME}-struct-placeholder",
            structs=structs,
            mode="struct",
        )


def get_app() -> PLF547App:
    """
    Get PLF547 benchmark app layout and callback registration.

    Returns
    -------
    PLF547App
        Benchmark layout and callback registration.
    """
    return PLF547App(
        name=BENCHMARK_NAME,
        title="PLF547",
        description=(
            "Performance in predicting protein-ligand fragment interaction energies "
            "for 547 protein fragment-ligand complexes."
        ),
        table_path=DATA_PATH / "plf547_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
            Div(id=f"{BENCHMARK_NAME}-struct-placeholder"),
        ],
    )


if __name__ == "__main__":
    # Create Dash app
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)

    # Construct layout and register callbacks
    plf547_app = get_app()
    full_app.layout = plf547_app.layout
    plf547_app.register_callbacks()

    # Run app
    full_app.run(port=8052, debug=True)
