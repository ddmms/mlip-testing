=============
Running tests
=============

This guide will break down how to run calculations, analysis, and the interactive
application.


Calculations
------------

Currently, all calculations should be launched using ``pytest``. This will help to
automatically discover and run each test, handle intermediate errors, and control
which tests are run based on our
`custom markers <https://docs.pytest.org/en/7.1.x/example/markers.html>`_.

ALl current tests can be launched using the
`run_calcs.sh <https://github.com/ddmms/ml-peg/blob/main/run_calcs.sh>`_ script:

.. code-block:: bash

    ./run_cals.sh


Individual tests or categories can be run using a similar command to the one in this
script, such as:

.. code-block:: bash

    pytest -v ml_peg/calcs/surfaces/*/calc* -s --run-slow


This will run all calculations in the surfaces category, including any marked as ``slow``.


Analysis
--------

As with calculations, analysis of results should also be launched using ``pytest``,
which can be done using the
`run_analysis.sh <https://github.com/ddmms/ml-peg/blob/main/run_analysis.sh>`_ script:

.. code-block:: bash

    ./run_analysis.sh


Individual tests or categories can be also analysed similarly. For example:

.. code-block:: bash

    pytest -v ml_peg/analysis/surfaces/*/analyse* -s


Will analyse the results of calculations in the surfaces category.


Application
-----------

Having run analysis, the app can now be launched by running the
`run_app.py <https://github.com/ddmms/ml-peg/blob/main/run_app.py>`_ Python script:

.. code-block:: bash

    python3 run_app.py

.. tip::

    ``uv run`` can be used in place of ``python3``, removing the need to activate a
    virtual environment.


By default, this will make the app visiable at http://localhost:8050.

.. tip::

    You can set the ``PORT`` environment variable to change the port used by Dash.


When launched, the app will attempt to automatically construct tables, figures, and
interactive features, based on any importable test apps defined in ``ml_peg/apps/``.

If any plots are unable to be loaded, a warning will be raised, and only the table will
be rendered for the test.

If a test's table is also unable to be loaded, the test will not be added to the app,
but the app builder should continue to attempt adding other tests.
