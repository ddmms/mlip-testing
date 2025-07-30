"""Base class to construct app layouts and register callbacks."""

from __future__ import annotations

from abc import ABC, abstractmethod

from dash.html import Div

from mlip_testing.app.utils.build import layout_builder


class BaseApp(ABC):
    """
    Abstract base class to construct app layouts and register callbacks.

    Parameters
    ----------
    title
        Title for application.
    table
        Dash table for application metrics.
    extra_components
        List of other Dash components to add to app.
    """

    def __init__(self, title, table, extra_components):
        """
        Initiaise class.

        Parameters
        ----------
        title
            Title for application.
        table
            Dash table for application metrics.
        extra_components
            List of other Dash components to add to app.
        """
        self.title = title
        self.table = table
        self.extra_components = extra_components
        self.layout = self.build_layout()

    def build_layout(self) -> Div:
        """
        Build app layout.

        Returns
        -------
        Div
            Div component with list all components for app.
        """
        # Define all components/placeholders
        return layout_builder(
            title=self.title, table=self.table, extra_components=self.extra_components
        )

    @abstractmethod
    def register_callbacks(self):
        """Register callbacks with app."""
        pass
