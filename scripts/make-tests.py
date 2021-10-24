#! /usr/bin/env python3 -B

from collections import OrderedDict as odict

def init(scene: odict):
  scene['asset'] = {}
  scene['cameras'] = []
  scene['environments'] = []
  scene['textures'] = []
  scene['shapes'] = []
  scene['subdivs'] = []
  scene['materials'] = []
  scene['instances'] = []

def find(array: list, name: str) -> int:
  for idx, element in enumerate(array):
    if element['name'] == name: return idx
  return -1

def asset(scene: odict, name: str) -> int :
  scene['asset'] = {
    "copyright": "Model by Fabio Pellacini from github.com/~xelatihy/yocto-gl",
    "generator": "Yocto/GL - https://github.com/xelatihy/yocto-gl",
    "version": "4.2"
  }
  return 0

def camera(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['cameras'], name) < 0:
    if name == "default":
      scene['cameras'] += [{
        "name": "default",
        "frame": [
          0.8151804208755493, -0.0, 0.579207181930542, 0.16660168766975403,
          0.9577393531799316, -0.23447643220424652, -0.5547295212745667,
          0.28763750195503235, 0.7807304263114929, -0.75, 0.4000000059604645,
          0.8999999761581421
        ],
        "aspect": 2.4000000953674316,
        "focus": 1.2168092727661133
      }]
    else:
      raise NameError(name)
  return find(scene['cameras'], name)

def texture(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['textures'], name) < 0:
    if name in ['floor', 'uvgrid', 'bumps-normal']:
      scene['textures'] += [{
        'name': name,
        'uri': 'textures/' + name + '.png'
      }]
    elif name in ['sky', 'sunsky']:
      scene['textures'] += [{
        'name': name,
        'uri': 'textures/' + name + '.hdr'
      }]
    else:
      raise NameError(name)
  return find(scene['textures'], name)

def environment(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['environments'], name) < 0:
    if name in ['sky', 'sunsky']:
      scene['environments'] += [{
        'name': name,
        'emission': [0.5, 0.5, 0.5],
        'emission_tex': texture(scene, name)
      }]
    else:
      raise NameError(name)
  return find(scene['environments'], name)

def shape(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['shapes'], name) < 0:
    if name in ['floor', 'sphere', 'bunny', 'arealight1', 'arealight2']:
      scene['shapes'] += [{
        'name': name,
        'uri': 'shapes/' + name + ".ply"
      }]
    else:
      raise NameError(name)
  return find(scene['shapes'], name)

def material(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['materials'], name) < 0:
    if name in ['floor']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'color': [1, 1, 1],
        'color_tex': texture(scene, name)
      }]
    elif name in ['arealight1', 'arealight2']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'emission': [20, 20, 20],
        'color': [0, 0, 0]
      }]
    elif name in ['uvgrid']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [1, 1, 1],
        'roughness': 0.2,
        'color_tex': texture(scene, name)
      }]
    elif name in ['jade']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [0.5, 0.5, 0.5],
        'roughness': 0,
        'scattering': [0.3, 0.6, 0.3]
      }]
    elif name in ['roughmetal']:
      scene['materials'] += [{
        'name': name,
        'type': 'reflective',
        'color': [0.66, 0.45, 0.34],
        'roughness': 0.2
      }]
    elif name in ['redglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [1.0, 0.5, 0.5],
        'roughness': 0
      }]
    elif name in ['bumped']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [0.5, 0.7, 0.5],
        'roughness': 0.2,
        'normal_tex': texture(scene, 'bumps-normal')
      }]
    else:
      raise NameError(name)
  return find(scene['materials'], name)

def instance(scene: odict, name: str, frame: list = []):
  if name == 'none': return -1
  if find(scene['instances'], name) < 0:
    shape_name, mat_name = name.split('-') if '-' in name else (name, name)
    scene['instances'] += [{
      'name': name,
      'frame': frame if frame else [1,0,0,0,1,0,0,0,1,0,0,0],
      'shape': shape(scene, shape_name),
      'material': material(scene, mat_name)
    }]
  return find(scene['instances'], name)

def instances(scene: odict, names: list, offset: list = [0, 0, 0], stride = 0.2):
  for idx, name in enumerate(names):
    origin = [ offset[0] + stride * (idx - len(names) // 2), offset[1], offset[2] ]
    frame = [1, 0, 0, 0, 1, 0, 0, 0, 1] + origin
    instance(scene, name, frame)

def arealights(scene: odict, name: str):
  if name == 'none': return -1
  instance(scene, 'arealight1', [
      0.8944271802902222, -0.0, 0.4472135901451111, 0.27562475204467773,
      0.7874992489814758, -0.5512495040893555, -0.3521803617477417,
      0.6163156628608704, 0.7043607234954834, -0.4000000059604645,
      0.800000011920929, 0.800000011920929])
  instance(scene, 'arealight2', [
      0.8944271802902222, 0.0, -0.4472135901451111, -0.27562475204467773,
      0.7874992489814758, -0.5512495040893555, 0.3521803617477417,
      0.6163156628608704, 0.7043607234954834, 0.4000000059604645,
      0.800000011920929, 0.800000011920929])

def test(scene: odict, name: str):
  init(scene)
  asset(scene, 'default')
  camera(scene, 'default')
  environment(scene, 'sky')
  arealights(scene, 'small')
  instance(scene, 'floor')
  if name in ['features1']:
    instances(scene, ['bunny-uvgrid', 'sphere-redglass', 'bunny-jade', 'sphere-bumped', 'bunny-roughmetal'])
  else:
    raise NameError(name)

def make_test(name: str) -> odict:
  scene = odict()
  test(scene, name)
  return scene

def make_scene(name: str, dirname: str = 'tests2'):
  import os, json, shutil
  scene = make_test(name)
  print(name)
  if os.path.exists(f'{dirname}/{name}'): shutil.rmtree(f'{dirname}/{name}')
  os.mkdir(f'{dirname}/{name}')
  os.mkdir(f'{dirname}/{name}/shapes')
  os.mkdir(f'{dirname}/{name}/textures')
  with open(f'{dirname}/{name}/{name}.json', 'w') as f:
    json.dump(scene, f, indent=2)
  for texture in scene['textures']:
    uri = texture['uri']
    shutil.copy(f'{dirname}/_assets/{uri}', f'{dirname}/{name}/{uri}')
  for shape in scene['shapes']:
    uri = shape['uri']
    shutil.copy(f'{dirname}/_assets/{uri}', f'{dirname}/{name}/{uri}')

if __name__ == '__main__':
  make_scene('features1')
