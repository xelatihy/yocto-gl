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
    'copyright': 'Model by Fabio Pellacini from github.com/~xelatihy/yocto-gl',
    'generator': 'Yocto/GL - https://github.com/xelatihy/yocto-gl',
    'version': '4.2'
  }
  return 0

def camera(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['cameras'], name) < 0:
    if name == 'default':
      scene['cameras'] += [{
        'name': 'default',
        'frame': [
          0.8151804208755493, -0.0, 0.579207181930542, 0.16660168766975403,
          0.9577393531799316, -0.23447643220424652, -0.5547295212745667,
          0.28763750195503235, 0.7807304263114929, -0.75, 0.4000000059604645,
          0.8999999761581421
        ],
        'aspect': 2.4000000953674316,
        'focus': 1.2168092727661133
      }]
    else:
      raise NameError(name)
  return find(scene['cameras'], name)

def texture(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['textures'], name) < 0:
    if name in ['floor', 'uvgrid', 'bumpsnormal', 'bumpsdisplacement']:
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
    if name in ['floor', 'sphere', 'cube', 'bunny', 'teapot', 'disk',
                'arealight1', 'arealight2', 'largearealight1', 'largearealight2', 
                'suzannesubdiv', 'displacedsubdiv', 'cubesubdiv',
                'hairball', 'hairballi',
                'flipcapuvsphere', 'uvcylinder']:
      scene['shapes'] += [{
        'name': name,
        'uri': 'shapes/' + name + '.ply'
      }]
    else:
      raise NameError(name)
  return find(scene['shapes'], name)

def subdiv(scene: odict, name: str) -> int:
  if name == 'none': return -1
  if find(scene['subdivs'], name) < 0:
    if name in ['suzannesubdiv']:
      scene['subdivs'] += [{
        'name': name,
        'shape': shape(scene, name),
        'uri': 'subdivs/' + name + '.obj',
        'subdivisions': 2,
        'catmullclark': True,
        'smooth': True
      }]
    elif name in ['displacedsubdiv']:
      scene['subdivs'] += [{
        'name': name,
        'shape': shape(scene, name),
        'uri': 'subdivs/' + name + '.obj',
        'catmullclark': True,
        'smooth': True,
        'displacement': 0.025,
        'displacement_tex': texture(scene, 'bumpsdisplacement')
      }]
    elif name in ['cubesubdiv']:
      scene['subdivs'] += [{
        'name': name,
        'shape': shape(scene, name),
        'uri': 'subdivs/' + name + '.obj',
        'subdivisions': 4,
        'catmullclark': True,
        'smooth': True
      }]
    else:
      raise NameError(name)
  return find(scene['subdivs'], name)

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
    elif name in ['largearealight1', 'largearealight2']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'emission': [10, 10, 10],
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
    elif name in ['matte']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'color': [0.7, 0.7, 0.7]
      }]
    elif name in ['plastic']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [0.5, 0.7, 0.5],
        'roughness': 0.2
      }]
    elif name in ['roughplastic']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [0.5, 0.7, 0.5],
        'roughness': 0.2
      }]
    elif name in ['sharpplastic']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [0.5, 0.5, 0.7],
        'roughness': 0
      }]
    elif name in ['glass']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [1.0, 1.0, 1.0],
        'roughness': 0
      }]
    elif name in ['jade']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [0.5, 0.5, 0.5],
        'roughness': 0,
        'scattering': [0.3, 0.6, 0.3]
      }]
    elif name in ['cloud']:
      scene['materials'] += [{
        'name': name,
        'type': 'volume',
        'color': [0.5, 0.5, 0.5],
        'roughness': 0,
        'scattering': [0.9, 0.9, 0.9]
      }]
    elif name in ['smoke']:
      scene['materials'] += [{
        'name': name,
        'type': 'volume',
        'color': [0.65, 0.65, 0.65],
        'roughness': 0,
        'scattering': [0.2, 0.2, 0.2]
      }]
    elif name in ['roughmetal']:
      scene['materials'] += [{
        'name': name,
        'type': 'reflective',
        'color': [0.66, 0.45, 0.34],
        'roughness': 0.2
      }]
    elif name in ['sharpmetal']:
      scene['materials'] += [{
        'name': name,
        'type': 'reflective',
        'color': [0.7, 0.7, 0.7],
        'roughness': 0
      }]
    elif name in ['redglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [1.0, 0.5, 0.5],
        'roughness': 0
      }]
    elif name in ['sharpglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [1.0, 1.0, 1.0],
        'roughness': 0
      }]
    elif name in ['roughglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'refractive',
        'color': [1.0, 0.7, 0.7],
        'roughness': 0.1
      }]
    elif name in ['sharpthinglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'transparent',
        'color': [1.0, 1.0, 1.0],
        'roughness': 0
      }]
    elif name in ['roughthinglass']:
      scene['materials'] += [{
        'name': name,
        'type': 'transparent',
        'color': [1.0, 0.7, 0.7],
        'roughness': 0.1
      }]
    elif name in ['notopaque']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'color': [0.7, 0.5, 0.5],
        'opacity': 0.2,
        'roughness': 0.1
      }]
    elif name in ['bumped']:
      scene['materials'] += [{
        'name': name,
        'type': 'glossy',
        'color': [0.5, 0.7, 0.5],
        'roughness': 0.2,
        'normal_tex': texture(scene, 'bumpsnormal')
      }]
    elif name in ['hair']:
      scene['materials'] += [{
        'name': name,
        'type': 'matte',
        'color': [0.7, 0.7, 0.7]
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
    if 'subdiv' in shape_name:
      subdiv(scene, shape_name)
  return find(scene['instances'], name)

def instances(scene: odict, names: list, offset: list = [0, 0, 0], stride = 0.2, interior: list = []):
  for idx, name in enumerate(names):
    origin = [ offset[0] + stride * (idx - len(names) // 2), offset[1], offset[2] ]
    frame = [1, 0, 0, 0, 1, 0, 0, 0, 1] + origin
    instance(scene, name, frame)
  for idx, name in enumerate(interior):
    if name == '': continue
    origin = [ offset[0] + stride * (idx - len(interior) // 2), offset[1], offset[2] ]
    frame = [1, 0, 0, 0, 1, 0, 0, 0, 1] + origin
    instance(scene, name, frame)

def arealights(scene: odict, name: str):
  if name == 'none': return -1
  if name == 'default':
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
  elif name in ['large']:
    instance(scene, 'largearealight1', [
        0.8944271802902222, -0.0, 0.4472135901451111, 0.2873478829860687,
        0.766261100769043, -0.5746957659721375, -0.3426823318004608,
        0.6425293684005737, 0.6853646636009216, -0.800000011920929,
        1.600000023841858, 1.600000023841858])
    instance(scene, 'largearealight2', [
        0.8944271802902222, 0.0, -0.4472135901451111, -0.2873478829860687,
        0.766261100769043, -0.5746957659721375, 0.3426823318004608,
        0.6425293684005737, 0.6853646636009216, 0.800000011920929,
        1.600000023841858, 1.600000023841858])
  else:
      raise NameError(name)

def test(scene: odict, name: str):
  init(scene)
  asset(scene, 'default')
  camera(scene, 'default')
  environment(scene, 'sky')
  if name in ['materials1', 'materials2', 'materials4']:
    arealights(scene, 'large')
  else:
    arealights(scene, 'default')
  instance(scene, 'floor')
  if name in ['features1']:
    instances(scene, ['bunny-uvgrid', 'sphere-redglass', 'bunny-jade', 'sphere-bumped', 'bunny-roughmetal'])
  elif name in ['features2']:
    instances(scene, ['sphere-uvgrid', 'suzannesubdiv-roughplastic', 'hairball-hair', 'displacedsubdiv-roughplastic', 'cube-uvgrid'],
      interior=['','','hairballi-hair','',''])
  elif name in ['materials1']:
    instances(scene, ['sphere-sharpplastic', 'sphere-roughplastic', 'sphere-matte', 'sphere-sharpmetal', 'sphere-roughmetal'])
  elif name in ['materials2']:
    instances(scene, ['sphere-sharpglass', 'sphere-roughglass', 'sphere-notopaque', 'sphere-sharpthinglass', 'sphere-roughthinglass'])
  elif name in ['materials4']:
    instances(scene, ['sphere-cloud', 'sphere-redglass', 'sphere-glass', 'sphere-jade', 'sphere-smoke'])
  elif name in ['shapes1']:
    instances(scene, ['sphere-uvgrid', 'flipcapuvsphere-uvgrid', 'disk-uvgrid', 'uvcylinder-uvgrid', 'cube-uvgrid'])
  elif name in ['shapes2']:
    instances(scene, ['cubesubdiv-uvgrid', 'suzannesubdiv-matte', 'displacedsubdiv-plastic', 'bunny-uvgrid', 'teapot-uvgrid'])
  else:
    raise NameError(name)

def make_test(name: str) -> odict:
  scene = odict()
  test(scene, name)
  return scene

def make_scene(name: str, dirname: str = 'tests'):
  import os, json, shutil
  scene = make_test(name)
  print(name)
  if os.path.exists(f'{dirname}/{name}'): shutil.rmtree(f'{dirname}/{name}')
  os.mkdir(f'{dirname}/{name}')
  if scene['shapes']: os.mkdir(f'{dirname}/{name}/shapes')
  if scene['textures']: os.mkdir(f'{dirname}/{name}/textures')
  if scene['subdivs']: os.mkdir(f'{dirname}/{name}/subdivs')
  with open(f'{dirname}/{name}/{name}.json', 'w') as f:
    json.dump(scene, f, indent=2)
  for texture in scene['textures']:
    uri = texture['uri']
    shutil.copy(f'{dirname}/_assets/{uri}', f'{dirname}/{name}/{uri}')
  for shape in scene['shapes']:
    uri = shape['uri']
    shutil.copy(f'{dirname}/_assets/{uri}', f'{dirname}/{name}/{uri}')
  for subdiv in scene['subdivs']:
    uri = subdiv['uri']
    shutil.copy(f'{dirname}/_assets/{uri}', f'{dirname}/{name}/{uri}')

if __name__ == '__main__':
  make_scene('features1')
  make_scene('features2')
  make_scene('materials1')
  make_scene('materials2')
  make_scene('materials4')
  make_scene('shapes1')
  make_scene('shapes2')
