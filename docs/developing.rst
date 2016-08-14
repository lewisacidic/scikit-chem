.. _developing:

Developing
==========

Development occurs on GitHub_. We gladly accept `pull requests`_ !

Development Requirements
------------------------

To start developing features for the package, you will need the core runtime
dependencies, shown in :ref:`installing <installing>`, in addition to the
below:

Testing
~~~~~~~

- py.test
- pytest-cov
- coverage

Linting
~~~~~~~

- pylint

Documentation
~~~~~~~~~~~~~

- sphinx >= 1.4
- sphinx_bootstrap_theme
- nbsphinx

These are all installable with pip.

Continuous Integration
----------------------

Pull requests and commits are automatically built and tested on Travis_.

Running the Tests
-----------------

Tests may be run locally through `py.test`_.  This can be invoked using either
``py.test`` or ``python setup.py test`` in the project root.  Command line
extensions are not tested by default - these can be tested also, by using the
appropriate flag, such as ``python setup.py test --with-chemaxon``.

Test Coverage
-------------

Test coverage is assessed using ``coverage``.  This is run locally as part of
the pytest command. It is set up to run as part of the CI, and can be viewed
on Scrutinzer_.  Test coverage has suffered as features were rapidly developed
in response to needs for the author's PhD, and will be improved once the PhD
is submitted!

Code Quality
------------

scikit-chem** conforms to pep8.  PyLint is used to assess code quality locally,
and can be run using ``pylint skchem`` from the root of the project.
Scrutinzer_ is also set up to run as part of the CI.  As with test coverage,
code quality has slipped due to time demands, and will be fixed once the PhD is
submitted!

Documentation
-------------

This documentation is built using Sphinx_, and Bootstrap using the Bootswatch
 Flatly theme.  The documentation is hosted on `Github Pages`_.  To build the
 html documentation locally, run ``make html``. To serve it, run
 ``make livehtml``.

.. _GitHub: https://github.com/richlewis42/scikit-chem
.. _pull requests: https://github.com/richlewis42/scikit-chem/pulls
.. _Github Pages: https://richlewis42.github.io/scikit-chem
.. _py.test: http://docs.pytest.org/en/latest/
.. _Scrutinzer: https://scrutinizer-ci.com/g/richlewis42/scikit-chem/
.. _Sphinx: http://www.sphinx-doc.org
.. _Travis: https://travis-ci.org/richlewis42/scikit-chem