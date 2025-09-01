"""Base class to construct app layouts and register callbacks."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from dash.development.base_component import Component
from dash.html import Div

from mlip_testing.app.utils.build_components import build_test_layout
from mlip_testing.app.utils.load import rebuild_table


class BaseApp(ABC):
    """
    Abstract base class to construct app layouts and register callbacks.

    Parameters
    ----------
    name
        Name for application tab.
    title
        Title for application.
    description
        Description of benchmark.
    table_path
        Path to json file containing Dash table data for application metrics.
    extra_components
        List of other Dash components to add to app.
    """

    def __init__(
        self,
        name: str,
        title: str,
        description: str,
        table_path: Path,
        extra_components: list[Component],
    ):
        """
        Initiaise class.

        Parameters
        ----------
        name
            Name for application test.
        title
            Title for benchmark.
        description
            Description of benchmark.
        table_path
            Path to json file containing Dash table data for application metrics.
        extra_components
            List of other Dash components to add to app.
        """
        self.name = name
        self.title = title
        self.description = description
        self.table_path = table_path
        self.extra_components = extra_components

        self.table_id = f"{self.name}-table"
        self.table = rebuild_table(self.table_path, id=self.table_id)
        self.layout = self.build_layout()

    def build_layout(self) -> Div:
        """
        Build layout for application.

        Returns
        -------
        Div
            Div component with list all components for app.
        """
        # Define all components/placeholders
        return build_test_layout(
            title=self.title,
            description=self.description,
            table=self.table,
            extra_components=self.extra_components,
        )

    @abstractmethod
    def register_callbacks(self):
        """Register callbacks with app."""
        pass
