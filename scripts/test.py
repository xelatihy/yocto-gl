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
    false = False # to cut and paste from json
    default_scene = {
        "name": "simple-ml",
        "cameras": [
            {
                "name": "default",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "!!proc": { "from": [-7.5, 4, 9], "to": [-0.75, 0.5, -0.5] }
                # "!!proc": { "from": [-3, 5, 10], "to": [0, 0, 0] }
            },
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
                "name": "matte",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
            },
            {
                "name": "plastic-sharp",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01
            },
            {
                "name": "plastic-rough",
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1
            },
            {
                "name": "metal-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0
            },
            {
                "name": "metal-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.1
            },
            {
                "name": "transparent",
                "diffuse": [ 0.7, 0.5, 0.5 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "opacity": 0.2
            },
            {
                "name": "glass-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": true
            },
            {
                "name": "glass-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.1,
                "refract": true
            },
            {
                "name": "thinglass-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": false
            },
            {
                "name": "thinglass-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.05,
                "refract": false
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
                "name": "teapot",
                "filename": "models/utah-teapot.obj"
            },
            {
                "name": "sphere",
                "filename": "models/test-sphere.obj",
                "!!proc": { "type": "sphere", "size": 1.5, "align_bottom": true }
            },
            {
                "name": "cube",
                "filename": "models/test-cube.obj",
                "!!proc": { "type": "cube_rounded", "size": [1.5, 1.5, 1.5], "align_bottom": true }
            },
            {
                "name": "disk",
                "filename": "models/test-disk.obj",
                "!!proc": { "type": "disk", "size": 1.5, "align_bottom": true }
            },
            {
                "name": "sphere-flipcap",
                "filename": "models/test-flipcap.obj",
                "!!proc": { "type": "uvsphere_flipcap", "size": 1.5, "align_bottom": true }
            },
            {
                "name": "cylinder",
                "filename": "models/test-cylinder.obj",
                "!!proc": { "type": "uvcylinder_rounded", "size": [1.5, 1.5, 1.5], "align_bottom": true }
            },
            {
                "name": "subdiv-cube",
                "filename": "models/test-subdiv-cube.obj",
                "subdivision_level": 4,
                "catmull_clark": true,
                "compute_normals": true,
                "preserve_facevarying": true,
                "!!proc": { "type": "cube_facevarying", "size": [1.5, 1.5, 1.5], "align_bottom": true }
            },
            {
                "name": "subdiv-monkey",
                "filename": "models/test-subdiv-monkey.obj",
                "subdivision_level": 2,
                "catmull_clark": true,
                "compute_normals": true,
                "!!proc": { "type": "suzanne", "size": 1.5, "align_bottom": true }
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
        import copy
        scene = copy.deepcopy(default_scene)
        scene['instances'] += copy.deepcopy(lights['instances'])
        scene['environments'] += copy.deepcopy(lights['environments'])
        posx = [ -4, -2, 0, 2, 4 ]
        for i in range(len(posx)):
            shape = shapes[i % len(shapes)]
            material = materials[i % len(materials)]
            scene['instances'] += [ {
                'name': f'{shape}_{material}',
                'shape': shape,
                'material': material,
                'frame': [ 1, 0, 0, 0, 1, 0, 0, 0, 1, posx[i], 0, 0 ]
            } ]
        with open(f'tests/{name}.json', 'wt') as f: json.dump(scene, f, indent=4)
    make_test('simple', ['bunny'], ['uvgrid'], area_lights)
    make_test('materials1', ['sphere'], ['plastic-sharp', 'plastic-rough', 'matte', 'metal-sharp', 'metal-rough'], area_lights)
    make_test('materials2', ['sphere'], ['glass-sharp', 'glass-rough', 'transparent', 'thinglass-sharp', 'thinglass-rough'], area_lights)
    make_test('materials3', ['sphere'], ['plastic-sharp-bumped', 'plastic-rough-bumped', 'matte-bumped', 'metal-sharp-bumped', 'metal-rough-bumped'], area_lights)
    make_test('shapes1', ['sphere', "sphere-flipcap", "disk", "cylinder", "cube"], ['uvgrid'], area_lights)
    make_test('shapes2', ['subdiv-cube', "subdiv-monkey", "teapot", "bunny", "subdiv-cube"], ['uvgrid', 'plastic-sharp'], area_lights)

cli()