#! /usr/bin/env python3 -B

import os, glob, json, sys
from typing import OrderedDict

def upgrade(filename, remove_names):
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
      nscene[groupname] = []
      for name, item in scene[groupname].items():
        nscene[groupname].append(OrderedDict())
        if not remove_names: nscene[groupname][-1]['name'] = name
        if groupname in ['textures', 'shapes']:
          nscene[groupname][-1]['uri'] = groupname + "/" + item
        else:
          nscene[groupname][-1].update(item)
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

parser = argparse.ArgumentParser(description='Upgrades scene.')
parser.add_argument('scenes', type=str, nargs='+')
parser.add_argument('--remove-names', action='store_true', default=False)
args = parser.parse_args()

for filename in args.scenes:
  print(filename)
  upgrade(filename, args.remove_names)
