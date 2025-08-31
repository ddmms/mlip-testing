"""Surface benchmark category node."""

from __future__ import annotations

from pathlib import Path

import zntrack

from mlip_testing.calcs.surfaces.OC157.calc_OC157 import OC157Benchmark
from mlip_testing.calcs.surfaces.S24.calc_S24 import S24Benchmark


class SurfaceCategory(zntrack.Node):
    """
    Category node for surface-related benchmarks.

    Aggregates OC157 and S24 benchmarks into a single category.
    """

    oc157_benchmarks: list[OC157Benchmark] = zntrack.deps()
    s24_benchmarks: list[S24Benchmark] = zntrack.deps()

    def run(self):
        """Organise the DAG (categories don't run calculations)."""
        pass


def build_surface_category_project(repro: bool = False) -> None:
    """
    Build surface category project with OC157 and S24 benchmarks.

    Parameters
    ----------
    repro
        Whether to call dvc repro -f after building.
    """
    import mlipx

    from mlip_testing.calcs.models.models import MODELS
    from mlip_testing.calcs.utils.utils import chdir

    project = mlipx.Project()

    # Create individual benchmarks
    oc157_benchmarks = []
    s24_benchmarks = []

    for model_name, model in MODELS.items():
        with project.group(model_name):
            oc157 = OC157Benchmark(
                model=model,
                model_name=model_name,
            )
            s24 = S24Benchmark(
                model=model,
                model_name=model_name,
            )
            oc157_benchmarks.append(oc157)
            s24_benchmarks.append(s24)

    # Create surface category
    with project.group("surface"):
        SurfaceCategory(
            oc157_benchmarks=oc157_benchmarks,
            s24_benchmarks=s24_benchmarks,
        )

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_surface_category():
    """Run surface category via pytest."""
    build_surface_category_project(repro=True)


if __name__ == "__main__":
    build_surface_category_project()
