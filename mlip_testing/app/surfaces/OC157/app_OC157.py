"""Run OC157 app."""

from __future__ import annotations

from pathlib import Path

from dash import Dash
from dash.html import Div
import numpy as np

from mlip_testing.app import APP_ROOT
from mlip_testing.app.base_app import BaseApp
from mlip_testing.app.utils.build_callbacks import (
    plot_from_table_column,
    struct_from_scatter,
)
from mlip_testing.app.utils.load import read_plot
from mlip_testing.calcs.models.models import MODELS

BENCHMARK_NAME = Path(__file__).name.removeprefix("app_").removesuffix(".py")
DATA_PATH = APP_ROOT / "data" / "surfaces" / "OC157"

SCATTER = read_plot(
    DATA_PATH / "figure_rel_energies.json", id=f"{BENCHMARK_NAME}-figure"
)
STRUCTS_DIR = DATA_PATH / list(MODELS.keys())[0]
# Assets dir will be parent directory
STRUCTS = list(
    np.repeat(
        [
            f"assets/surfaces/OC157/{list(MODELS.keys())[0]}/{i}.xyz"
            for i in range(len(list(STRUCTS_DIR.glob("*.xyz"))))
        ],
        3,
    )
)


class OC157App(BaseApp):
    """OC157 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={"MAE": SCATTER, "Ranking Error": SCATTER},
        )

        struct_from_scatter(
            scatter_id=f"{BENCHMARK_NAME}-figure",
            struct_id=f"{BENCHMARK_NAME}-struct-placeholder",
            structs=STRUCTS,
            mode="traj",
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
        name=BENCHMARK_NAME,
        title="OC157",
        description=(
            "Performance in predicting relative energies between 3 structures for 157 "
            "molecule-surface combinations."
        ),
        table_path=DATA_PATH / "oc157_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
            Div(id=f"{BENCHMARK_NAME}-struct-placeholder"),
        ],
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
