.. _installing:

Installation and Getting Started
================================

``scikit-chem`` is easy to install and configure.  Detailed instructions are
listed below.  The quickest way to get everything installed is by
`using conda`_.

.. todo::

    fish

Installation
------------

``scikit-chem`` is tested on Python 2.7 and 3.5. It depends on rdkit_, most of
the core `Scientific Python Stack`_, as well as several smaller pure *Python*
libraries.

Dependencies
~~~~~~~~~~~~

The full list of dependencies is:

- rdkit_
- numpy_
- scipy_
- matplotlib_
- scikit-learn_
- pandas_
- h5py_
- fuel_
- ipywidgets_
- progressbar2_

The package and these dependencies are available through two different *Python*
package managers, conda_ and pip_.  It is recommended to use conda_.

Using conda_
~~~~~~~~~~~~

conda_ is a cross-platform, Python-agnostic package and environment manager
developed by `Continuum Analytics`_.  It provides packages as prebuilt binary
files, allowing for straightforward installation of Python packages, even those
with complex C/C++ extensions.  It is installed as part of the Anaconda_
Scientific Python distribution, or as the lightweight miniconda_.

The package and all dependencies for ``scikit-chem`` are available through the
defaults_ or richlewis_ conda_ channel.  To install:

.. code:: bash

  conda install -c richlewis scikit-chem


This will install ``scikit-chem`` with all its dependencies from the
`author's anaconda repository`_ as ``conda`` packages.

.. attention::

    For Windows, you will need to install a dependency, fuel_, separately.
    This will be made available via conda_ in the future.


Using pip_
~~~~~~~~~~

pip_ is the standard Python package manager.  The package is available via
PyPI_, although the dependencies may require compilation or at worst may not
work at all.

.. code:: bash

    pip install scikit-chem

This will install ``scikit-chem`` with all available dependencies as regular
``pip`` controlled packages.

.. attention::

    A key dependency, rdkit_, is not installable using pip_, and will need to
    be installed by other means, such as conda_, apt-get_ on Linux or Homebrew_
    on Mac, or compiled from source (not recommended!!).

Configuration
-------------

scikit-chem
~~~~~~~~~~~

Currently, scikit-chem cannot be configured in a config file.  This feature is
planned to be added in future releases.  To request this feature as a priority,
please mention it in the appropriate
`github issue <https://github.com/richlewis42/scikit-chem/issues/7>`_

Fuel
~~~~

To use the :ref:`data <tutorial/data>` functionality, you will need to set up
fuel_.  This involves configuring the .fuelrc.  An example .fuelrc might be as
follows:

.. code-block:: yaml

    data_path: ~/datasets
    extra_downloaders:
    - skchem.data.downloaders
    extra_converters:
    - skchem.data.converters

This adds the location for fuel datasets, and adds the ``scikit-chem`` data
downloaders and converters to the fuel command line tools.

.. _conda: http://conda.pydata.org/docs/index.html
.. _Continuum Analytics: https://www.continuum.io
.. _Scientific Python Stack: https://www.scipy.org/stackspec.html
.. _Anaconda: https://www.continuum.io/downloads
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _numpy: http://numpy.org
.. _scipy: http://scipy.org
.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _scikit-learn: https://scikit-learn.org
.. _progressbar2: http://progressbar-2.readthedocs.io/en/latest/
.. _ipywidgets: https://ipywidgets.readthedocs.io/en/latest/
.. _h5py: http://h5py.org
.. _`author's anaconda repository`: https://conda.anaconda.org/richlewis
.. _richlewis: https://conda.anaconda.org/richlewis
.. _defaults: https://conda.anaconda.org/defaults
.. _fuel: https://github.com/mila-udem/fuel
.. _PyPI: https://pypi.python.org/pypi
.. _pip: https://pypi.python.org/pypi/pip
.. _apt-get: http://linux.die.net/man/8/apt-get
.. _Homebrew: http://brew.sh
.. _rdkit: http://www.rdkit.org