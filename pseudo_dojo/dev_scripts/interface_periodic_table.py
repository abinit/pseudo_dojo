#!/usr/bin/env python
import sys
import os

import matplotlib.pyplot as plt
from ptplotter.plotter import ElementDataPlotter

extension = 'psp8'

table = {}

for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(extension):
             element =  os.path.join(root).split('/')[-1].split('-')[-1].capitalize()
             print(element)
             try:
                 element_data = table[element]
             except KeyError:
                 element_data = {'n': 0}
             element_data['n'] += 1
             table[element] = element_data

print(table)

plotter = ElementDataPlotter(data=table)

def np(elt):
    """Number of data sets"""
    try:
        return elt['n']
    except KeyError:
        return 0

plotter.ptable(functions=[np], cmaps='RdYlGn', color='Black')
plt.show()

