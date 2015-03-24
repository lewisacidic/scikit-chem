#! /bin/bash
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variables defined
# in the .travis.yml in the top level folder of the project.

set -e

# update travis environment
sudo apt-get update

# retrieve miniconda distribution appropriate for travis' python version
if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
else
  wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

# install miniconda
bash miniconda.sh -b -p $HOME/miniconda

# add it to path
export PATH="$HOME/miniconda/bin:$PATH"

# configure conda to not ask for confirmation and not to change the prompt
conda config --set always_yes yes --set changeps1 no

# update conda (shouldn't be necessary as the latest is retrieved, but can't hurt)
conda update -q conda

# Useful for debugging any issues with conda
conda info -a

# Add required channels for dependencies
conda config --add channels 'http://conda.binstar.org/rdkit'

# Create the virtual environment
conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION 

# Activate the virtual environment
source activate test-environment

# Install the package and test dependencies
conda install --file $TRAVIS_BUILD_DIR/requirements.txt
conda install --file $TRAVIS_BUILD_DIR/test_requirements.txt

# Install the package
python $TRAVIS_BUILD_DIR/setup.py install