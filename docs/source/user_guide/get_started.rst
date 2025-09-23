===============
Getting started
===============

Dependencies
------------

All required and optional dependencies can be found in `pyproject.toml <https://github.com/ddmms/ml-peg/blob/main/pyproject.toml>`_.


Installation
------------

The latest stable release of ``ML-PEG``, including its dependencies, will soon be installable from PyPI by running:

.. code-block:: bash

    python3 -m pip install ml-peg


To get all the latest changes, ``ML-PEG`` can also be installed from GitHub:

.. code-block:: bash

    python3 -m pip install git+https://github.com/ddmms/ml-peg.git


Running the application
-----------------------

To build the container using the ``Dockerfile`` provided, run:

.. code-block:: bash

    docker build . -t ml-peg-app


Once built, you can mount your current directory and start the app by running:

.. code-block:: bash

    docker run --volume .:/app  --publish 8050:8050 benchmark-app


Alternatively, you can use the ``compose.yml`` file provided, via Docker Compose:

.. code-block:: bash

    docker compose up -d


The app should now be accessible at http://localhost:8050.
