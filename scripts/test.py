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
    default_scene = {
        "name": "simple-ml",
        "cameras": [
            {
                "name": "cam",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.0576,
                "film_height": 0.024,
                "!!proc": {
                    "from": [ -3, 5, 10 ],
                    "to": [ 0, 0, 0 ]
                },
            }
        ],
        "textures": [
            {
                "name": "floor",
                "filename": "textures/floor.png",
                "gamma": 2.2,
                "!!proc": { "type": "grid", "border": true }
            },
            {
                "name": "uvgrid",
                "filename": "textures/uvgrid.png",
                "gamma": 2.2,
                "!!proc": { "type": "uvgrid" }
            }
        ],
        "materials": [
            {
                "name": "floor",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_texture": "floor"
            },
            {
                "name": "uvgrid",
                "diffuse": [ 1, 1, 1 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1,
                "diffuse_texture": "uvgrid"
            },
            {
                "name": "arealight1",
                "emission": [ 20, 20, 20 ]
            },
            {
                "name": "arealight2",
                "emission": [ 20, 20, 20 ]
            }
        ],
        "shapes": [
            {
                "name": "floor",
                "filename": "models/floor.ply",
                "!!proc": { "type": "floor" }
            },
            {
                "name": "bunny",
                "filename": "models/stanford-bunny.obj"
            },
            {
                "name": "test-sphere",
                "filename": "models/obj2.ply",
                "!!proc": { "type": "sphere", "size": 1.5, "align_bottom": true }
            },
            {
                "name": "test-cube",
                "filename": "models/obj3.ply",
                "!!proc": { "type": "cube_rounded", "size": [1.5, 1.5, 1.5], "align_bottom": true }
            },
            {
                "name": "arealight1",
                "!!proc": { "type": "quad", "size": [ 4, 4 ] }
            },
            {
                "name": "arealight2",
                "!!proc": { "type": "quad", "size": [ 4, 4 ] }
            }
        ],
        "instances": [
            {
                "name": "floor",
                "shape": "floor",
                "material": "floor"
            }
        ],
        "environments": []
    }
    area_lights = {
        "instances": [
            {
                "name": "arealight1",
                "shape": "arealight1",
                "material": "arealight1",
                "!!proc": { "from": [ -4, 8, 8 ], "to": [ 0, 1, 0 ] }
            },
            {
                "name": "arealight2",
                "shape": "arealight2",
                "material": "arealight2",
                "!!proc": { "from": [ 4, 8, 8 ], "to": [ 0, 1, 0 ] }
            }
        ],
        "environments": []
    }
    def make_test(name, shapes, materials, lights):
        scene = default_scene.copy()
        scene['instances'] += lights['instances']
        scene['environments'] += lights['environments']
        posx = [ -4, -2, 0, 2, 4, -4, -2, 0, 2, 4 ]
        posy = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
        posz = [ 2, 2, 2, 2, 2, -2, -2, -2, -2, -2 ]
        for i in range(10):
            shape = shapes[i % len(shapes)]
            material = materials[i % len(materials)]
            scene['instances'] += [ {
                'name': f'{shape}_{material}',
                'shape': shape,
                'material': material,
                'frame': [ 1, 0, 0, 0, 1, 0, 0, 0, 1, posx[i], posy[i], posz[i] ]
            } ]
        with open(f'tests/{name}.json', 'wt') as f: json.dump(scene, f, indent=4)
    make_test('simple', ['bunny'], ['uvgrid'], area_lights)

cli()