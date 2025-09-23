=============
Running Tests
=============

This guide will break down how to run calculations, analysis, and the application.


Calculations
++++++++++++

Currently, all calculations should be launched using ``pytest``. This will help to
automatically discover and run each test, as well as providing decorators such as
``fixture``, ``mark.parametrize``, and custom ``mark`` options.

ALl current tests can be launched using
`run_calcs.sh <https://github.com/ddmms/ML-PEG/blob/main/run_calcs.sh>`_.

Individual tests or categories can be run using a similar command to the one in this
script, such as:

.. code-block:: bash

    pytest -v ml_peg/calcs/surfaces/*/calc* -s --run-slow


This will run all calculations in the surfaces category, including any marked as ``slow``.


Analysis
++++++++

As with calculations, analysis of results should also be launched using ``pytest``,
which can be done using the
`run_analysis.sh <https://github.com/ddmms/ML-PEG/blob/main/run_analysis.sh>`_ script.

Individual tests or categories can be also analysed similarly. For example:

.. code-block:: bash

    pytest -v ml_peg/analysis/surfaces/*/analyse* -s


Will analyse the results of calculations in the surfaces category.


Application
+++++++++++

Having run analysis, the app can now be launched by running
`run_app.py <https://github.com/ddmms/ML-PEG/blob/main/run_app.py>`_:

.. code-block:: bash

    python3 run_app.py


By default, this will make the app visiable at http://localhost:8050.

.. tip::

    You can set the ``PORT`` environment variable to change the port used by Dash.
