#!/usr/bin/env python
import sys
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser

path = sys.argv[1]

parser = OncvOutputParser(path)
parser.scan()

print(parser.core)
print(parser.valence)
print(parser.rc_min)

#1s 2s 2p 3s 3p 3d 4s 4p 4d
#5s 5p 4f 5d 6s
#1.4

