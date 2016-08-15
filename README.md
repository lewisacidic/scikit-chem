[![Travis](https://img.shields.io/travis/richlewis42/scikit-chem.svg?style=flat-square)](https://travis-ci.org/richlewis42/scikit-chem)
[![Coverage](https://img.shields.io/scrutinizer/coverage/g/richlewis42/scikit-chem.svg?style=flat-square)](https://scrutinizer-ci.com/g/richlewis42/scikit-chem)
[![Code Health](https://img.shields.io/scrutinizer/g/richlewis42/scikit-chem.svg?style=flat-square)](https://scrutinizer-ci.com/g/richlewis42/scikit-chem/)
[![Release](https://img.shields.io/pypi/v/scikit-chem.svg?style=flat-square)](https://github.com/richlewis42/scikit-chem/releases)
[![PyPI](https://img.shields.io/pypi/dm/scikit-chem.svg?style=flat-square)](https://pypi.python.org/pypi/scikit-chem)
[![Conda](https://anaconda.org/richlewis/scikit-chem/badges/installer/conda.svg)](https://anaconda.org/richlewis/scikit-chem)
[![DOI](https://zenodo.org/badge/4513/richlewis42/scikit-chem.svg?style=flat-square)](https://zenodo.org/badge/latestdoi/4513/richlewis42/scikit-chem)

# scikit-chem

**scikit-chem** is a high level cheminformatics library, built on the excellent [RDKit](https://github.com/rdkit/rdkit), aimed to integrating with the Scientific Python Stack.
scikit-chem is a high level cheminformatics library built on rdkit that aims to integrate with the Scientific Python Stack by promoting interoperativity with libraries such as pandas and scikit-learn, and emulating similar patterns and APIs as found in those libraries.

Some notable features include:

- Pythonic core API
- Consistent, declarative interfaces for many cheminformatics tasks, including:
  - Reading file formats
  - Chemical standardization
  - Conformer generation
  - Filtering
  - Feature calculation
  - Pipelining
- A simple interface for chemical datasets
- Structure visualization
- Interactivity in Jupyter Notebooks

scikit-chem should be thought of as a simple complement to the excellent rdkit - scikit-chem objects are subclasses of rdkit objects, and as such, the two libraries can usually be used together easily when the advanced functionality of rdkit is required.
