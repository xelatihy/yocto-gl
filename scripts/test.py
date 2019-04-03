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
        os.system('rm tests/_results/*.png')
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
                "!!proc": { "from": [-0.75, 0.4, 0.9], "to": [-0.075, 0.05, -0.05] }
            },
            {
                "name": "front",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "!!proc": { "from": [0, 0.575, 0.14], "to": [0, 0.05, 0] }
            },
            {
                "name": "back",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "!!proc": { "from": [0, 0.575, -0.14], "to": [0, 0.05, 0] }
            },
            {
                "name": "perspective-sharp",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "!!proc": { "from": [-0.75, 0.4, 0.9], "to": [-0.075, 0.05, -0.05] }
            },
            {
                "name": "perspective-dof",
                "focal_length": 0.05,
                "lens_aperture": 0.25,
                "film_width": 0.036,
                "film_height": 0.015,
                "!!proc": { "from": [-0.75, 0.4, 0.9], "to": [-0.075, 0.05, -0.05] }
            },
            {
                "name": "orthographic-sharp",
                "focal_length": 0.005,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": true,
                "!!proc": { "from": [-0.75, 0.4, 0.9], "to": [-0.075, 0.05, -0.05] }
            },
            {
                "name": "orthographic-dof",
                "focal_length": 0.005,
                "lens_aperture": 0.2,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": true,
                "!!proc": { "from": [-0.75, 0.4, 0.9], "to": [-0.075, 0.05, -0.05] }
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
            },
            {
                "name": "bumps",
                "filename": "textures/bumps.png",
                "!!proc": { "type": "bump", "bump_scale": 0.05 }
            },
            {
                "name": "bumps-normal",
                "filename": "textures/bumps-normal.png",
                "!!proc": { "type": "bump", "bump_to_normal": true, "bump_scale": 0.05 }
            },
            {
                "name": "fbm-displacement",
                "filename": "textures/fbm-displacement.png",
                "height_scale": 0.025,
                "!!proc": { "type": "fbm", "scale": 10 }
            },
            {
                "name": "sky",
                "filename": "textures/sky.hdr",
                "!!proc": { "type": "sky" }
            },
            {
                "name": "sunsky",
                "filename": "textures/sky.hdr",
                "!!proc": { "type": "sky", "has_sun": true, "width": 2048, "height": 1024 }
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
                "name": "matte-displaced",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "displacement_texture": "fbm-displacement"
            },
            {
                "name": "plastic-sharp-bumped",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01,
                "normal_texture": "bumps-normal"
            },
            {
                "name": "metal-sharp-bumped",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0,
                "normal_texture": "bumps-normal"
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
                "name": "hair",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
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
                "filename": "models/test-floor.ply",
                "!!proc": { "type": "floor", "size": [ 4, 4 ] }
            },
            {
                "name": "bunny",
                "filename": "models/test-bunny.obj"
            },
            {
                "name": "teapot",
                "filename": "models/test-teapot.obj"
            },
            {
                "name": "sphere",
                "filename": "models/test-sphere.obj",
                "!!proc": { "type": "sphere", "size": 0.15, "align_bottom": true }
            },
            {
                "name": "cube",
                "filename": "models/test-cube.obj",
                "!!proc": { "type": "box_rounded", "size": [0.15, 0.15, 0.15], "rounded": 0.3, "align_bottom": true }
            },
            {
                "name": "disk",
                "filename": "models/test-disk.obj",
                "!!proc": { "type": "disk", "size": 0.15, "align_bottom": true }
            },
            {
                "name": "sphere-flipcap",
                "filename": "models/test-flipcap.obj",
                "!!proc": { "type": "uvsphere_flipcap", "size": 0.15, "align_bottom": true }
            },
            {
                "name": "cylinder",
                "filename": "models/test-cylinder.obj",
                "!!proc": { "type": "uvcylinder_rounded", "size": [0.15, 0.15, 0.15], "align_bottom": true }
            },
            {
                "name": "sphere-displaced",
                "filename": "models/test-sphere-displaced.obj",
                "preserve_facevarying": false,
                "!!proc": { "type": "sphere", "size": 0.15, "align_bottom": true }
            },
            {
                "name": "subdiv-cube",
                "filename": "models/test-subdiv-cube.obj",
                "subdivision_level": 4,
                "catmull_clark": true,
                "compute_normals": true,
                "preserve_facevarying": true,
                "!!proc": { "type": "cube_facevarying", "size": [0.15, 0.15, 0.15], "align_bottom": true }
            },
            {
                "name": "subdiv-monkey",
                "filename": "models/test-subdiv-monkey.obj",
                "subdivision_level": 2,
                "catmull_clark": true,
                "compute_normals": true,
                "!!proc": { "type": "suzanne", "size": 0.15, "align_bottom": true }
            },
            {
                "name": "hairball1",
                "filename": "models/test-hairball1.ply",
                "!!proc": { "type": "hairball", "size": 0.15, "noise": [ 0.03, 100 ] }
            },
            {
                "name": "hairball2",
                "filename": "models/test-hairball2.ply",
                "!!proc": { "type": "hairball", "size": 0.15 }
            },
            {
                "name": "hairball3",
                "filename": "models/test-hairball3.ply",
                "!!proc": { "type": "hairball",  "size": 0.15, "clump": [ 0.5, 128 ] }
            },
            {
                "name": "hairballi",
                "filename": "models/test-hairballi.ply",
                "!!proc": { "type": "hairball_interior", "size": 0.15 }
            },
            {
                "name": "arealight1",
                "!!proc": { "type": "quad", "size": [ 0.4, 0.4 ] }
            },
            {
                "name": "arealight2",
                "!!proc": { "type": "quad", "size": [ 0.4, 0.4 ] }
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
    mixed_lights = {
        "instances": [
            {
                "name": "arealight1",
                "shape": "arealight1",
                "material": "arealight1",
                "!!proc": { "from": [ -0.4, 0.8, 0.8 ], "to": [ 0, 0.1, 0 ] }
            },
            {
                "name": "arealight2",
                "shape": "arealight2",
                "material": "arealight2",
                "!!proc": { "from": [ 0.4, 0.8, 0.8 ], "to": [ 0, 0.1, 0 ] }
            }
        ],
        "environments": [
            {
                "name": 'sky',
                "emission": [2, 2, 2],
                "emission_texture": "sky"
            }
        ]
    }
    sunsky_lights = {
        "instances": [],
        "environments": [
            {
                "name": 'sunsky',
                "emission": [1, 1, 1],
                "emission_texture": "sunsky"
            }
        ]
    }
    def make_test(name, shapes, materials, lights, xoffsets=[ -0.4, -0.2, 0, 0.2, 0.4 ], yoffsets=[0,0,0,0,0], zoffsets=[0,0,0,0,0], xscales=[1,1,1,1,1], yscales=[1,1,1,1,1], zscales=[1,1,1,1,1]):
        import copy
        scene = copy.deepcopy(default_scene)
        scene['instances'] += copy.deepcopy(lights['instances'])
        scene['environments'] += copy.deepcopy(lights['environments'])
        num = max(len(xoffsets), len(shapes), len(materials))
        for i in range(num):
            shape = shapes[i % len(shapes)]
            material = materials[i % len(materials)]
            xoffset = xoffsets[i % len(xoffsets)]
            yoffset = yoffsets[i % len(yoffsets)]
            zoffset = zoffsets[i % len(zoffsets)]
            xscale = xscales[i % len(xscales)]
            yscale = yscales[i % len(yscales)]
            zscale = zscales[i % len(zscales)]
            if not shape: continue
            if not material: continue
            scene['instances'] += [ {
                'name': f'{shape}_{material}',
                'shape': shape,
                'material': material,
                'frame': [ xscale, 0, 0, 0, yscale, 0, 0, 0, zscale, xoffset, yoffset, zoffset ]
            } ]
        if True: # pruning
            old_materials = scene['materials']
            scene['materials'] = []
            for material in old_materials:
                used = False
                for instance in scene['instances']:
                    if instance['material'] == material['name']: used = True
                if used: scene['materials'] += [material] 
            old_shapes = scene['shapes']
            scene['shapes'] = []
            for shape in old_shapes:
                used = False
                for instance in scene['instances']:
                    if instance['shape'] == shape['name']: used = True
                if used: scene['shapes'] += [shape] 
            old_textures = scene['textures']
            scene['textures'] = []
            for texture in old_textures:
                used = False
                for material in scene['materials']:
                    if 'emission_texture' in material and material['emission_texture'] == texture['name']: used = True
                    if 'diffuse_texture' in material and material['diffuse_texture'] == texture['name']: used = True
                    if 'normal_texture' in material and material['normal_texture'] == texture['name']: used = True
                    if 'displacement_texture' in material and material['displacement_texture'] == texture['name']: used = True
                for environment in scene['environments']:
                    if environment['emission_texture'] == texture['name']: used = True
                if used: scene['textures'] += [texture] 
        with open(f'tests/{name}.json', 'wt') as f: json.dump(scene, f, indent=4)
    make_test('features1', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ['uvgrid', 'plastic-sharp', 'metal-rough', 'plastic-rough', 'metal-sharp'], mixed_lights)
    make_test('materials1', ['sphere'], ['plastic-sharp', 'plastic-rough', 'matte', 'metal-sharp', 'metal-rough'], mixed_lights)
    make_test('materials2', ['sphere'], ['glass-sharp', 'glass-rough', 'transparent', 'thinglass-sharp', 'thinglass-rough'], mixed_lights)
    make_test('materials3', ['sphere', 'sphere', 'sphere-displaced', 'sphere', 'sphere'], ['plastic-sharp-bumped', 'plastic-sharp-bumped', 'matte-displaced', 'metal-sharp-bumped', 'metal-sharp-bumped'], mixed_lights)
    make_test('shapes1', ['sphere', "sphere-flipcap", "disk", "cylinder", "cube"], ['uvgrid'], mixed_lights)
    make_test('shapes2', ['subdiv-cube', "subdiv-monkey", "teapot", "bunny", "subdiv-cube"], ['uvgrid', 'plastic-sharp'], mixed_lights)
    make_test('shapes3', ['sphere', "hairball1", "hairball2", "hairball3", "sphere", "", "hairballi", "hairballi", "hairballi", ""], ['matte', 'hair', 'hair', 'hair', 'matte'], mixed_lights,
        yoffsets=[ 0, 0.075, 0.075, 0.075, 0 ], xscales=[ 0.5, 1, 1, 1, 0.5 ])
    make_test('arealights1', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ['uvgrid', 'plastic-sharp', 'metal-rough', 'plastic-rough', 'metal-sharp'], area_lights)
    make_test('environments1', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ['uvgrid', 'plastic-sharp', 'metal-rough', 'plastic-rough', 'metal-sharp'], sunsky_lights)

cli()