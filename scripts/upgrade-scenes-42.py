#! /usr/bin/env python3 -B

import os, glob, json, sys
from typing import OrderedDict

def upgrade(filename):
  with open(filename) as f:
    scene = json.load(f, object_pairs_hook=OrderedDict)
    if scene['asset']['version'] == "4.2": return
    if scene['asset']['version'] != "4.1": return
    shape_map = { name: id for id, name in enumerate(scene['shapes'].keys()) } if 'shapes' in scene else {}
    texture_map = { name: id for id, name in enumerate(scene['textures'].keys()) } if 'textures' in scene else {}
    material_map = { name: id for id, name in enumerate(scene['materials'].keys()) } if 'materials' in scene else {}
    nscene = OrderedDict()
    nscene['asset'] = scene['asset']
    nscene['asset']['version'] = "4.2"
    for groupname in ['cameras', 'environments', 'textures', 'materials', 'shapes', 'subdivs', 'instances']:
      if groupname not in scene: continue
      if groupname in ['textures', 'shapes']:
        nscene[groupname] = [ { 'uri': groupname + "/" + item } for _, item in scene[groupname].items() ]
      else:
        nscene[groupname] = [ item for _, item in scene[groupname].items() ]
    if 'subdivs' in nscene:
      for subdiv in nscene['subdivs']:
        subdiv['uri'] = 'subdivs/' + subdiv['datafile']
        del subdiv['datafile']
        subdiv['shape'] = shape_map[subdiv['shape']]
    for material in nscene['materials']:
      for key in material:
        if key.endswith('_tex'): material[key] = texture_map[material[key]]
    if 'environments' in nscene:
      for environment in nscene['environments']:
        for key in environment:
          if key.endswith('_tex'): environment[key] = texture_map[environment[key]]
    for instance in nscene['instances']:
      instance['shape'] = shape_map[instance['shape']]
      instance['material'] = material_map[instance['material']]
  with open(filename, 'w') as f:
    json.dump(nscene, f, indent=2)

for filename in sys.argv[1:]:
  print(filename)
  upgrade(filename)
