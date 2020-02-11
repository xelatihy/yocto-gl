import glob, json, os

for filename in glob.glob('tests/*.json'):
  print(filename)
  with open(filename) as f: scene = json.load(f)
  scene['materials'] = []
  for shape in scene['shapes']:
    material = {}
    material['name'] = shape['name'].replace('shapes/', 'materials/')
    keys = [ key for key in shape if key not in ['shape','instances','name','frame','lookat'] ]
    for key in keys:
      material[key] = shape[key]
      del shape[key]
    shape['material'] = material['name']
    scene['materials'] += [ material ]
  with open(filename, 'w') as f: json.dump(scene, f, indent=2)
