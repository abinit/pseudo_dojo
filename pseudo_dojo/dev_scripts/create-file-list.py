#!/usr/bin/env python

import os
import glob
import json

jsons = glob.glob('*.json')
tgzs = glob.glob('pseudos/*.tgz')
pseudos = glob.glob('pseudos/*/*')

files = jsons+tgzs+pseudos

with open('files.json', 'w') as f:
    json.dump(files, fp=f)

exit()
