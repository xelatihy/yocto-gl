#! /usr/bin/env python3 -B

import click, glob, os, json

@click.group()
def cli():
    pass

@cli.command()
def run():
    os.system('mkdir -p build && mkdir -p build/release && cd build/release && cmake ../.. -GNinja -DYOCTO_EMBREE=ON')
    os.system('rm tests/_output/*.png; rm  tests/_difference/*.png')
    os.system('cd build/release && ctest -j 4 --output-on-failure')

@cli.command()
def clean():
    os.system('rm tests/_output/*.png; rm tests/_difference/*.png')

@cli.command()
@click.option('--clean/--no-clean', default=False)
def update(clean=False):
    if clean:
        os.system('rm tests/_result/*.png')
    os.system('cp tests/_output/*.png tests/_results')

@cli.command()
def format():
    from collections import OrderedDict
    for filename in sorted(glob.glob('tests/*.json')):
        with open(filename, 'r') as f: js = json.load(f, object_pairs_hook=OrderedDict)
        with open(filename, 'w') as f: json.dump(js, f, indent=4)

@cli.command()
def make_tests():
    true = True   # to cut and paste from json
    flase = False # to cut and paste from json
    default_cameras = {
        "cameras": [
            {
                "name": "cam",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "!!proc": {
                    "from": [ -3, 5, 10 ],
                    "to": [ 0, 0, 0 ]
                },
                "film_width": 0.0576,
                "film_height": 0.024
            }
        ]
    }
    default_floor = {
        "textures": [
            {
                "name": "floor",
                "filename": "textures/floor.png",
                "gamma": 2.2,
                "!!proc": { "type": "grid", "border": true }
            }
        ],
        "materials": [
            {
                "name": "floor",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_texture": "floor"
            }
        ],
        "shapes": [
            {
                "name": "floor",
                "filename": "models/floor.ply",
                "!!proc": { "type": "floor" }
            }
        ],
        "instances": [
            {
                "name": "floor",
                "shape": "floor",
                "material": "floor"
            },
        ]
    }
    area_lights = {
        "materials": [
            {
                "name": "light1",
                "emission": [ 20, 20, 20 ]
            },
            {
                "name": "light2",
                "emission": [ 20, 20, 20 ]
            }
        ],
        "shapes": [
            {
                "name": "light1",
                "!!proc": { "type": "quad", "size": [ 4, 4 ] }
            },
            {
                "name": "light2",
                "!!proc": { "type": "quad", "size": [ 4, 4 ] }
            }
        ],
        "instances": [
            {
                "name": "light1",
                "shape": "light1",
                "material": "light1",
                "!!proc": { "from": [ -4, 8, 8 ], "to": [ 0, 1, 0 ] }
            },
            {
                "name": "light2",
                "shape": "light2",
                "material": "light2",
                "!!proc": { "from": [ 4, 8, 8 ], "to": [ 0, 1, 0 ] }
            }
        ]
    }
    default_materials = {
        "textures": [
            {
                "name": "uvgrid",
                "filename": "textures/uvgrid.png",
                "gamma": 2.2,
                "!!proc": { "type": "uvgrid" }
            }
        ],
        "materials": [
            {
                "name": "uvgrid",
                "diffuse": [ 1, 1, 1 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1,
                "diffuse_texture": "uvgrid"
            },
        ]
    }
    default_shapes = {
        "shapes": [
            {
                "name": "test-sphere",
                "filename": "models/test-sphere.ply",
                "!!proc": { "type": "sphere", "size": 1.5, "align_bottom": true }
            }
        ]
    }
    def make_test(name, scenes, make_instances=True):
        js = {
            'name': name,
            'cameras': [],
            'textures': [],
            'materials': [],
            'shapes': [],
            'instances': [],
            'environments': [],
        }
        for scene in scenes:
            for key in scene:
                js[key] += scene[key]
        if make_instances:
            posx = [ -4, -2, 0, 2, 4, -4, -2, 0, 2, 4 ]
            posy = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
            posz = [ 2, 2, 2, 2, 2, -2, -2, -2, -2, -2 ]
            materials = scenes[-2]['materials']
            shapes = scenes[-1]['shapes']
            for i in range(10):
                js['instances'] += [ {
                    'name': f'obj{i}',
                    'shape': shapes[i % len(shapes)]['name'],
                    'material': materials[i % len(materials)]['name'],
                    'frame': [ 1, 0, 0, 0, 1, 0, 0, 0, 1, posx[i], posy[i], posz[i] ]
                } ]
        with open(f'tests/{name}.json', 'wt') as f: json.dump(js, f, indent=2)
    make_test('simple', [default_cameras, default_floor, area_lights, default_materials, default_shapes])

cli()