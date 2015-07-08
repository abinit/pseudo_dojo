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
# pseudos = multi.change_icmod3(fcfact_list=(3,), rcfact_list=(1.3, 1.35))
pseudos = multi.change_icmod3()
#pseudos = multi.change_icmod3(fcfact_list=(3, 4, 5), rcfact_list=(1.2, 1.25, 1.6, 1.65, 1.7, 1.75))

commands = []
for i, p in enumerate(pseudos):
#    print("[%i] %s" % (i, p.filepath))
    dir = os.path.split(os.path.split(p.filepath)[0])[-1]
    psp8 = os.path.split(p.filepath)[-1] 
    log = 'log ' #os.path.join(dir, "log")
    cmd = "nohup dojorun.py %s --trials=df --paral-kgb=0 > %s &" % (psp8, psp8 + '.log')
#    cmd = "cd %s\ncp ../*.yml .\nnohup dojorun.py %s --trials=df --paral-kgb=0 > %s &\ncd .." % (dir, psp8, log)
    commands.append(cmd)

for cmd in commands:
    print(cmd)
