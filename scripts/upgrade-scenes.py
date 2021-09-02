#! /usr/bin/env python3 -B

import os, glob, json, sys
from typing import OrderedDict

def upgrade(filename):
  with open(filename) as f:
    scene = json.load(f, object_pairs_hook=OrderedDict)
    scene['asset']['version'] = "4.1"
    scene['shapes'] = OrderedDict()
    for _, instance in scene['instances'].items():
      name = instance['shape']
      scene['shapes'][name] = name + ".ply"
    scene['textures'] = OrderedDict()
    for _, material in scene['materials'].items():
      if 'type' in material:
        if material['type'] == 'metallic':
          material['type'] = 'reflective'
      else:
        material['type'] = 'matte'
      for key in ['color_tex', 'normal_tex', 'emission_tex']:
        if key not in material: continue
        name = material[key]
        scene['textures'][name] = name + ".png"
    if 'environments' in scene:
      for _, environment in scene['environments'].items():
        for key in ['emission_tex']:
          if key not in environment: continue
          name = environment[key]
          scene['textures'][name] = name + ".hdr"
    # for _, subdiv in scene['subdivs'].items():
    #   subdiv['filename'] = name + ".obj"
    if not scene['textures']: del scene['textures']
  with open(filename, 'w') as f:
    json.dump(scene, f, indent=2)

for filename in sys.argv[1:]:
  print(filename)
  upgrade(filename)
