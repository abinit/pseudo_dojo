#!/usr/bin/env python

__author__ = 'setten'

import os
import sys
import collections
import numpy as np
import pprint

from pymatgen.util.convergence import determine_convergence
from abipy.flowtk.netcdf import NetcdfReader
from abipy.flowtk.pseudos import Pseudo


class DeltaFactorData(object):
    """
    class for delta factor data analisys
    """

    def __init__(self, ps_name):
        """

        """
        self.param = {}
        self.df_data = {}
        self.results = {}
        self.etotal_data = {}
        self.df_extra = np.inf
        self.ps_name = ps_name
        self.pseudo = Pseudo.from_file(ps_name+'.psp8')

    def read(self):
        """

        """
        tree = os.walk(self.ps_name+'_df_run_full')
        for dirs in tree:
            file_name = os.path.join(dirs[0], 'deltadata.txt')
            if os.path.isfile(file_name):
                f = open(file_name, 'r')
                lines = f.readlines()
                try:
                    df = float(lines[0].split()[3])
                except ValueError:
                    print('warning', lines[0].split()[3])
                    df = abs(complex(lines[0].split()[3]))
                #print lines[0], df
                #print dirs[0]
                location = os.path.split(dirs[0])
                #print location
                base = os.path.join(*location[0:len(location)-1])
                out_file = os.path.join(base, 't3', 'outdata', 'out_GSR.nc')
                #print out_file
                out = NetcdfReader(out_file)
                if not isinstance(out.read_value('ecut'), collections.Iterable):
                    ecut = out.read_value('ecut')
                    etotal = out.read_value('etotal')
                elif not isinstance(out.read_value('ecut')[0], collections.Iterable):
                    ecut = out.read_value('ecut')[0]
                    etotal = out.read_value('etotal')[0]
                else:
                    raise Exception
                self.etotal_data.update({ecut: etotal})
                self.df_data.update({ecut: df})
        self.pseudo.dojo_report['delta_factor'] = self.df_data
        self.pseudo.dojo_report['total_energy'] = self.etotal_data

    def print_data(self):
        """

        """
        pprint.pprint(self.df_data)

    def print_results(self, stream=sys.stdout):
        """

        """
        table = [["Convergense criterium", "ecut"]]

        print('extrapolated Delta Factor: ', my_df_data.df_extra)
        for x in sorted(self.results.keys()):
            print(x,  self.results[x])

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
            test_res = determine_convergence(xs, ys, 'df'+str(abs(tol)), tol=tol, verbose=False, mode='extra_noise')
            self.results.update({abs(tol): test_res[1]})
            self.df_extra = test_res[4]
        pprint.pprint(self.results)
        self.pseudo.dojo_report.update({'hints': {'high': self.results[1.0], 'medium': self.results[3.0],
                                                  'low': self.results[10], 'based_on': 'delta_factor'}})


if __name__ == "__main__":
    name = sys.argv[1]
    my_df_data = DeltaFactorData(name)
    my_df_data.read()
    my_df_data.test_convergence()
    pprint.pprint(my_df_data.pseudo.dojo_report)
    my_df_data.pseudo.write_dojo_report()
    my_df_data.print_results()

