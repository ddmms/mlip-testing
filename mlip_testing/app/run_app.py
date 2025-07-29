"""Run main application."""

from __future__ import annotations

from importlib import import_module
import os
from pathlib import Path

from dash import Dash
from dash.html import Div

from mlip_testing import app

DATA_PATH = Path(__file__).parent / "data"


def build_tabs(full_app: Dash) -> list[Div]:
    """
    Get layout and register callbacks for all tab applications.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.

    Returns
    -------
    list[Div]
        Layouts for all tabs.
    """
    # Find Python files e.g. app_OC157.py in mlip_tesing.app module.
    tabs = Path(app.__file__).parent.glob("*/app*.py")
    layout = []

    # Build all layouts, and register all callbacks to main app.
    for tab in tabs:
        tab_name = (tab.parent).name
        tab_app = import_module(f"mlip_testing.app.{tab_name}.app_{tab_name}")
        layout.append(tab_app.build_layout())
        tab_app.register_callbacks(full_app)

    return layout


def build_full_app(full_app: Dash):
    """
    Build full app layout and register callbacks.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.
    """
    full_app.layout = build_tabs(full_app)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8050))

    full_app = Dash(__name__, assets_folder=DATA_PATH)
    build_full_app(full_app)

    print(f"Starting Dash app on port {port}...")
    full_app.run(host="0.0.0.0", port=port)
