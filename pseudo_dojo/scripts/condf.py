__author__ = 'setten'

import os
import sys
import numpy as np
import collections

from pymatgen.io.gwwrapper.convergence import test_conv
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.util.string_utils import pprint_table

class DeltaFactorData(object):
    """
    class for delta factor data analisys
    """

    def __init__(self):
        """

        """
        self.param = {}
        self.df_data = {}
        self.results = {}
        self.df_extra = np.inf

    def read(self):
        """

        """
        tree = os.walk('df_run')
        for dirs in tree:
            file_name = os.path.join(dirs[0], 'deltadata.txt')
            if os.path.isfile(file_name):
                f = open(file_name, 'r')
                lines = f.readlines()
                df = float(lines[0].split()[3])
                #print lines[0], df
                #print dirs[0]
                location = os.path.split(dirs[0])
                #print location
                base = os.path.join(*location[0:len(location)-1])
                out_file = os.path.join(base, 't0', 'outdata', 'out_OUT.nc')
                #print out_file
                out = NetcdfReader(out_file)
                if not isinstance(out.read_value('ecut'), collections.Iterable):
                    ecut = out.read_value('ecut')
                elif not isinstance(out.read_value('ecut')[0], collections.Iterable):
                    ecut = out.read_value('ecut')[0]
                else:
                    raise Exception
                self.df_data.update({ecut: df})

    def print_data(self):
        """

        """
        print self.df_data

    def print_results(self, stream=sys.stdout):
        """

        """
        table = [["Convergense criterium", "ecut"]]

        print sorted(self.results)

        for x in sorted(self.results.keys()):

            table.append([x, self.results[x]])

        print 'extrapolated Delta Factor: ', my_df_data.df_extra
        print table
        pprint_table(table, out=stream)

    def test_convergence(self):
        xs = sorted(self.df_data.keys())
        #print xs
        ys = []
        for x in xs:
            ys.append(self.df_data[x])
        #print ys
        for tol in [-0.1, -0.01, -0.001, -0.0001]:
            test_res = test_conv(xs, ys, 'df', tol=tol, verbose=False)
            self.results.update({abs(tol): test_res[1]})
            self.df_extra = test_res[4]

if __name__ == "__main__":
    my_df_data = DeltaFactorData()
    my_df_data.read()
    my_df_data.test_convergence()
    my_df_data.print_results()

