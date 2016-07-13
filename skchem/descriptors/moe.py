#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.moe

Module specifying moe descriptors.
"""

class MOEDescriptorCalculator(object):


    def __init__(self):
        pass

    def transform(self, obj):
        if isinstance(obj, core.Mol):
            return self._transform_series(pd.Series(obj)).iloc[0]
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_series(self, series):

        with tempfile.NamedTemporaryFile(suffix='.sdf') as in_file, tempfile.NamedTemporaryFile() as out_file:
            # write mols to file
            write_sdf(series, in_file.name)
            args = ['mddesc', in_file.name, '-o', out_file.name] + self.index

            LOGGER.info('Running: ' + ' '.join(args))

            # call command line
            subprocess.call(args)
            try:
                finished = pd.read_table(out_file.name).set_index('id')
            except Exception:
                finished = None
        finished.index = series.index
        return finished


