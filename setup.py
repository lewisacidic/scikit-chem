from setuptools import setup

setup(name='scikit-chem', 
    version='0.0.1', 
    description='Python Scientific Stack integrated Cheminformatics library, using RDKit',
    url='http://github.com/richlewis42/scikit-chem',
    author='Richard Lewis',
    author_email='rl403@cam.ac.uk',
    license='BSD',
    packages=['skchem', 'skchem.ext'],
    package_data = {'skchem.ext': ['data/PIDGIN_models.pkl.gz']},
    install_requires=[
        #'rdkit',  ## currently rdkit does not have a valid pip installable package, so this dependency must be met by the installer
        'pandas',
        'numpy',
        'scikit-learn'
        ],
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe=False
    )
