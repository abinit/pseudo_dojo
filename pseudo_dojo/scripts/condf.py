__author__ = 'setten'

import os
import sys
import numpy as np
import collections

from pymatgen.io.gwwrapper.convergence import test_conv
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.util.string_utils import pprint_table
from pymatgen.io.abinitio.pseudos import Pseudo

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
        self.pseudo = Pseudo.from_file('totest')

    def read(self):
        """

        """
        tree = os.walk('df_run_full')
        for dirs in tree:
            file_name = os.path.join(dirs[0], 'deltadata.txt')
            if os.path.isfile(file_name):
                f = open(file_name, 'r')
                lines = f.readlines()
                try:
                    df = float(lines[0].split()[3])
                except ValueError:
                    print 'warning', lines[0].split()[3]
                    df = abs(complex(lines[0].split()[3]))
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
        self.pseudo.dojo_report.update({'delta_factor': self.df_data})

    def print_data(self):
        """

        """
        print self.df_data

    def print_results(self, stream=sys.stdout):
        """

        """
        table = [["Convergense criterium", "ecut"]]

        print 'extrapolated Delta Factor: ', my_df_data.df_extra
        for x in sorted(self.results.keys()):
            print x,  self.results[x]

        #    table.append([x, self.results[x]])

        #print table
        #pprint_table(table, out=stream)

    def test_convergence(self):
        xs = sorted(self.df_data.keys())
        #print xs
        ys = []
        for x in xs:
            ys.append(self.df_data[x])
        #print ys
        for tol in [-10, -3.0, -1.0, -0.3, -0.1]:
            test_res = test_conv(xs, ys, 'df'+str(abs(tol)), tol=tol, verbose=False, mode='extra_noise')
            self.results.update({abs(tol): test_res[1]})
            self.df_extra = test_res[4]
        print self.results
        self.pseudo.dojo_report.update({'hints': {'high': self.results[1.0], 'medium': self.results[3.0],
                                                  'low': self.results[10], 'based_on': 'delta_factor'}})


if __name__ == "__main__":
    my_df_data = DeltaFactorData()
    my_df_data.read()
    my_df_data.test_convergence()
    print my_df_data.pseudo.dojo_report
    my_df_data.pseudo.write_dojo_report()
    my_df_data.print_results()

