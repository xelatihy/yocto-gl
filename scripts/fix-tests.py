#! /usr/bin/env python3 -B

import os, glob, json

for filename in glob.glob('tests/_version40/*/*.json'):
  output = filename.replace('_version40/','')
  os.system(f'./bin/yscene convert {filename} --output {output}')
