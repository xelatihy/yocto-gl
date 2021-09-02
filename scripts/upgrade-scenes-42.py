#! /usr/bin/env python3 -B

import json, argparse
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
        nitem = OrderedDict()
        nscene[groupname].append(nitem)
        if not remove_names: nitem['name'] = name
        if groupname in ['textures', 'shapes']:
          nitem['uri'] = groupname + "/" + item
        else:
          nitem.update(item)
        for key, value in nitem.items():
          if key.endswith('_tex'): nitem[key] = texture_map[value]
          if key == 'shape': nitem[key] = shape_map[value]
          if key == 'material': nitem[key] = material_map[value]
          if key == 'frame' and isinstance(value[0], list):
            nitem[key] = value[0] + value[1] + value[2] + value[3]
        if 'datafile' in nitem:
          nitem['uri'] = item['datafile']
          del item['datafile']
  with open(filename, 'w') as f:
    json.dump(nscene, f, indent=2)

parser = argparse.ArgumentParser(description='Upgrades scene.')
parser.add_argument('scenes', type=str, nargs='+')
parser.add_argument('--remove-names', action='store_true', default=False)
args = parser.parse_args()

for filename in args.scenes:
  print(filename)
  upgrade(filename, args.remove_names)
