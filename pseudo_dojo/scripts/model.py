#!/usr/bin/env python
import sys
import os
from pseudo_dojo.ppcodes.ppgen import OncvMultiGenerator

#usage: model.py input_file

# Build the object that allows us to change the parameters in the template
path = sys.argv[1]
multi = OncvMultiGenerator(path)

# Change the parameters of the model core charge, generate pseudos
# and save results in new directories.
# HEre we generate 2 pseudos, default args are: fcfact_list=(3, 4, 5), rcfact_list=(1.3, 1.35, 1.4, 1.45, 1.5, 1.55))
pseudos = multi.change_icmod3(fcfact_list=(3,), rcfact_list=(1.3, 1.35))

commands = []
for i, p in enumerate(pseudos):
    print("[%i] %s" % (i, p.filepath))
    log = os.path.join(os.path.dirname(p.filepath), "log")
    cmd = "nohup dojorun.py %s --trials=df --paral-kgb=0 > %s &" % (p.filepath, log)
    commands.append(cmd)

for cmd in commands:
    print(cmd)
