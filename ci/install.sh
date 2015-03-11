#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e

sudo apt-get update

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
else
  wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda

# Useful for debugging any issues with conda
conda info -a

# Add required channels for dependencies
conda config --add channels 'http://conda.binstar.org/rdkit'

# Create the virtual environment with appropriate requirements
conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION --file $TRAVIS_BUILD_DIR/requirements.txt

# Activate the virtual environment
source activate test-environment

# Install the package
python $TRAVIS_BUILD_DIR/setup.py install