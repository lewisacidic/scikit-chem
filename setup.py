#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from setuptools.command.test import test as TestCommand
from setuptools import setup, find_packages
import sys

DISTNAME = 'scikit-chem'
DESCRIPTION = 'A set of python modules for cheminformatics'
with open('README.md') as f:
    LONG_DESCRIPTION = f.read()
MAINTAINER = 'Richard Lewis'
MAINTAINER_EMAIL = 'rl403@cam.ac.uk'
URL = 'http://github.com/richlewis42/scikit-chem'
LICENSE = 'new BSD'
DOWNLOAD_URL = 'https://github.com/richlewis42/scikit-chem/archive/master.zip'
CLASSIFIERS = [ 
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Education',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering :: Chemistry'
]
MAJOR = 0
MINOR = 0
MICRO = 6
VERSION = '{major}.{minor}.{micro}'.format(major=MAJOR, minor=MINOR, micro=MICRO)

with open('requirements.txt') as f:
    REQUIREMENTS = [l.strip() for l in f]
    REQUIREMENTS.remove('rdkit')

with open('test_requirements.txt') as f:
    TEST_REQUIREMENTS = [l.strip() for l in f]


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

def setup_package():

    metadata = dict(name=DISTNAME,
                    maintainer=MAINTAINER,
                    maintainer_email=MAINTAINER_EMAIL,
                    description=DESCRIPTION,
                    license=LICENSE,
                    url=URL,
                    version=VERSION,
                    download_url=DOWNLOAD_URL,
                    long_description=LONG_DESCRIPTION,
                    classifiers=CLASSIFIERS
    )

    setup(
        packages=find_packages(),
        cmdclass={'test': PyTest},
        package_data = {
            'skchem.target_prediction': ['data/PIDGIN_models.pkl.gz'],
            'skchem.resource': ['atom_data.csv'],
            'skchem.standardizers': ['default_config.xml']},
        include_package_data=True,
        install_requires=REQUIREMENTS,
        tests_require=TEST_REQUIREMENTS,
        zip_safe=False,
        **metadata
        )


if __name__ == "__main__":
    setup_package()
