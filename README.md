# ml-peg

[![Build Status][ci-badge]][ci-link]
[![Coverage Status][cov-badge]][cov-link]
[![Docs status][docs-badge]][docs-link]
[![License][license-badge]][license-link]
[![DOI][doi-badge]][doi-link]

Testing framework for machine learnt interatomic potentials.

The interactive analysis suite is currently hosted at: http://mlip-testing.stfc.ac.uk:8050

> [!NOTE]
> Migration in progress! The live benchmarks are currently run and analysed using
> [mlipx](https://github.com/basf/mlipx) nodes defined in this repository:
> https://github.com/joehart2001/mlipx.
>
> New benchmarks are expected to be added following the format defined in this
> repository, and work is ongoing to migrate all existing benchmarks to this format.

## Contents
- [Getting started](#getting-started)
- [Features](#features)
- [Development](#development)
- [License](#license)

## Getting started

### Dependencies

All required and optional dependencies can be found in [pyproject.toml](pyproject.toml).


### Installation

The latest stable release of `ml-peg`, including its dependencies, will be installable from PyPI by running:

```
python3 -m pip install ml-peg
```

To get all the latest changes, `ml-peg` can be installed from GitHub:

```
python3 -m pip install git+https://github.com/ddmms/ml-peg.git
```

## Features

Coming soon!


## Development

Please ensure you have consulted our
[contribution guidelines](https://github.com/ddmms/ml-peg/blob/main/contributing.md)
and
[coding style](https://github.com/ddmms/ml-peg/blob/main/coding_style.md)
before proceeding.

We recommend installing `uv` for dependency management when developing for `ml-peg`:

1. Install [uv](https://docs.astral.sh/uv/getting-started/installation)
2. Install `ml-peg` with dependencies in a virtual environment:

```shell
git clone https://github.com/ddmms/ml-peg
cd ml-peg
uv sync # Create a virtual environment and install dependencies
source .venv/bin/activate
pre-commit install  # Install pre-commit hooks
pytest -v  # Discover and run all tests
```

Please refer to the [online documentation](https://ddmms.github.io/ml-peg/developer_guide/index.html)
for information about contributing new benchmarks and models.

## License

[GNU General Public License version 3](LICENSE)

[ci-badge]: https://github.com/ddmms/ml-peg/actions/workflows/ci.yml/badge.svg?branch=main
[ci-link]: https://github.com/ddmms/ml-peg/actions
[cov-badge]: https://coveralls.io/repos/github/ddmms/ml-peg/badge.svg?branch=main
[cov-link]: https://coveralls.io/github/ddmms/ml-peg?branch=main
[docs-badge]: https://github.com/ddmms/ml-peg/actions/workflows/docs.yml/badge.svg
[docs-link]: https://ddmms.github.io/ml-peg/
[license-badge]: https://img.shields.io/badge/License-GPLv3-blue.svg
[license-link]: https://opensource.org/license/gpl-3-0
[doi-link]: https://doi.org/10.5281/zenodo.16904445
[doi-badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.16904445.svg
