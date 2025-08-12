"""Run main application."""

from __future__ import annotations

import os
from pathlib import Path

from dash import Dash

from mlip_testing.app.build_app import build_full_app

DATA_PATH = Path(__file__).parent / "data"


def run_app():
    """Set port and run Dash application."""
    port = int(os.environ.get("PORT", 8050))

    full_app = Dash(__name__, assets_folder=DATA_PATH)
    build_full_app(full_app)

    print(f"Starting Dash app on port {port}...")
    full_app.run(host="0.0.0.0", port=port, debug=True)


if __name__ == "__main__":
    run_app()
