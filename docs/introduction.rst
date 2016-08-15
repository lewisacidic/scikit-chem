.. _introduction:

An introduction to scikit-chem
==============================

**scikit-chem** is a high level cheminformatics library built on rdkit_ that
aims to integrate with the `Scientific Python Stack`_ by promoting
interoperativity with libraries such as pandas_ and scikit-learn_, and
emulating similar patterns and APIs as found in those libraries.

Some notable features include:

- *Pythonic* core API
- Consistent, declarative interfaces for many cheminformatics tasks, including:
     - Reading file formats
     - Chemical standardization
     - Conformer generation
     - Filtering
     - Feature calculation
     - Pipelining
- A simple interface for chemical datasets
- Structure visualization
- Interactivity in `Jupyter Notebooks`_

**scikit-chem** should be thought of as a simple complement to the excellent
rdkit_ - **scikit-chem** objects are subclasses of rdkit_ objects, and as such,
the two libraries can usually be used together easily when the advanced
functionality of rdkit_ is required.

.. _rdkit: http://www.rdkit.org
.. _pandas: http://pandas.pydata.org
.. _Scientific Python Stack: http://www.scipy.org
.. _scikit-learn: https://scikit-learn.org
.. _Jupyter Notebooks: http://jupyter.org
