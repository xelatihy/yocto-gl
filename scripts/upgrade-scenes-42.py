#! /usr/bin/env python3 -B

import os, glob, json, sys
from typing import OrderedDict

def upgrade(filename, named_refs):
  with open(filename) as f:
    scene = json.load(f, object_pairs_hook=OrderedDict)
    if scene['asset']['version'] == "4.2": return
    if scene['asset']['version'] != "4.1": return
    nscene = OrderedDict()
    nscene['asset'] = scene['asset']
    nscene['asset']['version'] = "4.2"
    for groupname in ['cameras', 'environments', 'textures', 'materials', 'shapes', 'subdivs', 'instances']:
      if groupname in ['textures', 'shapes']:
        nscene[groupname] = [ { 'uri': groupname + "/" + item } for item in scene[groupname] ]
      else:
        nscene[groupname] = [ item for item in scene[groupname] ]
    if 'subdivs' in nscene:
      for subdiv in nscene['subdivs']:
        subdiv['uri'] = 'subdivs/' + subdiv['datafile']
        del subdiv['datafile']
  with open(filename, 'w') as f:
    json.dump(nscene, f, indent=2)

for filename in sys.argv[1:]:
  print(filename)
  upgrade(filename)
