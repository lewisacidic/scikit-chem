
import os
import pandas as pd


def resource(*args):

    """ passes a file path for a data resource specified """

    return os.path.join(os.path.dirname(__file__), *args)

PERIODIC_TABLE = pd.read_csv(resource('atom_data.csv'), index_col=0)
ORGANIC = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
