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

def distance(eye, center):
    def length(a):
        from math import sqrt
        return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    return length([eye[0] - center[0], eye[1] - center[1], eye[2] - center[2]])

def lookat(eye, center, up, flipped=False):
    def length(a):
        from math import sqrt
        return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    def normalize(a):
        l = length(a)
        return [ a[0] / l, a[1] / l, a[2] / l ]
    def cross(a, b):
        return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

    w = normalize([eye[0] - center[0], eye[1] - center[1], eye[2] - center[2]])
    u = normalize(cross(up, w))
    v = normalize(cross(w, u))
    if flipped:
        w = [-w[0], -w[1], -w[2] ]
        u = [-u[0], -u[1], -u[2] ]
    return [ u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2], eye[0], eye[1], eye[2] ]

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
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "front",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, 0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, 1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "name": "back",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, -0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, -1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "name": "perspective-sharp",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "perspective-dof",
                "focal_length": 0.05,
                "lens_aperture": 0.025,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "orthographic-sharp",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": true,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "orthographic-dof",
                "focal_length": 0.05,
                "lens_aperture": 0.02,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": true,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
        ],
        "textures": [
            {
                "name": "test-floor",
                "filename": "textures/test-floor.png.ypreset"
            },
            {
                "name": "test-uvgrid",
                "filename": "textures/test-uvgrid.png.ypreset"
            },
            {
                "name": "test-bump",
                "filename": "textures/test-bump.png.ypreset"
            },
            {
                "name": "test-bump-normal",
                "filename": "textures/test-bump-normal.png.ypreset"
            },
            {
                "name": "test-fbm-displacement",
                "filename": "textures/test-fbm-displacement.png.ypreset"
            },
            {
                "name": "test-sky",
                "filename": "textures/test-sky.hdr.ypreset"
            },
            {
                "name": "test-sunsky",
                "filename": "textures/test-sunsky.hdr.ypreset"
            }
        ],
        "materials": [
            {
                "name": "test-floor",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_texture": "test-floor"
            },
            {
                "name": "test-uvgrid",
                "diffuse": [ 1, 1, 1 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1,
                "diffuse_texture": "test-uvgrid"
            },
            {
                "name": "test-matte",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
            },
            {
                "name": "test-plastic-sharp",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01
            },
            {
                "name": "test-plastic-rough",
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1
            },
            {
                "name": "test-metal-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0
            },
            {
                "name": "test-metal-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.1
            },
            {
                "name": "test-matte-displaced",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "displacement_texture": "test-fbm-displacement",
                "displacement_scale": 0.025
            },
            {
                "name": "test-plastic-sharp-bumped",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01,
                "normal_texture": "test-bump-normal"
            },
            {
                "name": "test-metal-sharp-bumped",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0,
                "normal_texture": "test-bump-normal"
            },
            {
                "name": "test-transparent",
                "diffuse": [ 0.7, 0.5, 0.5 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "opacity": 0.2
            },
            {
                "name": "test-glass-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": true
            },
            {
                "name": "test-glass-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.1,
                "refract": true
            },
            {
                "name": "test-thinglass-sharp",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": false
            },
            {
                "name": "test-thinglass-rough",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.05,
                "refract": false
            },
            {
                "name": "test-hair",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
            },
            {
                "name": "test-arealight1",
                "emission": [ 20, 20, 20 ]
            },
            {
                "name": "test-arealight2",
                "emission": [ 20, 20, 20 ]
            }
        ],
        "shapes": [
            {
                "name": "test-floor",
                "filename": "shapes/test-floor.ply.ypreset"
            },
            {
                "name": "test-bunny",
                "filename": "shapes/test-bunny.obj"
            },
            {
                "name": "test-teapot",
                "filename": "shapes/test-teapot.obj"
            },
            {
                "name": "test-sphere",
                "filename": "shapes/test-sphere.obj.ypreset"
            },
            {
                "name": "test-cube",
                "filename": "shapes/test-cube.obj.ypreset"
            },
            {
                "name": "test-disk",
                "filename": "shapes/test-disk.obj.ypreset"
            },
            {
                "name": "test-uvsphere-flipcap",
                "filename": "shapes/test-uvsphere-flipcap.obj.ypreset"
            },
            {
                "name": "test-uvcylinder",
                "filename": "shapes/test-uvcylinder.obj.ypreset"
            },
            {
                "name": "test-sphere-displaced",
                "filename": "shapes/test-sphere-displaced.obj.ypreset",
                "preserve_facevarying": false
            },
            {
                "name": "test-cube-subdiv",
                "filename": "shapes/test-cube-subdiv.obj.ypreset",
                "subdivision_level": 4,
                "catmull_clark": true,
                "compute_normals": true,
                "preserve_facevarying": true
            },
            {
                "name": "test-suzanne-subdiv",
                "filename": "shapes/test-suzanne-subdiv.obj.ypreset",
                "subdivision_level": 2,
                "catmull_clark": true,
                "compute_normals": true
            },
            {
                "name": "test-hairball1",
                "filename": "shapes/test-hairball1.ply.ypreset"
            },
            {
                "name": "test-hairball2",
                "filename": "shapes/test-hairball2.ply.ypreset"
            },
            {
                "name": "test-hairball3",
                "filename": "shapes/test-hairball3.ply.ypreset"
            },
            {
                "name": "test-hairball-interior",
                "filename": "shapes/test-hairball-interior.ply.ypreset"
            },
            {
                "name": "test-arealight1",
                "filename": "shapes/test-arealight1.ply.ypreset"
            },
            {
                "name": "test-arealight2",
                "filename": "shapes/test-arealight2.ply.ypreset"
            }
        ],
        "instances": [
            {
                "name": "test-floor",
                "shape": "test-floor",
                "material": "test-floor"
            }
        ],
        "environments": []
    }
    area_lights = {
        "instances": [
            {
                "name": "test-arealight1",
                "shape": "test-arealight1",
                "material": "test-arealight1",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "name": "test-arealight2",
                "shape": "test-arealight2",
                "material": "test-arealight2",
                "frame": lookat([ 0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            }
        ],
        "environments": []
    }
    mixed_lights = {
        "instances": [
            {
                "name": "test-arealight1",
                "shape": "test-arealight1",
                "material": "test-arealight1",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "name": "test-arealight2",
                "shape": "test-arealight2",
                "material": "test-arealight2",
                "frame": lookat([ 0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            }
        ],
        "environments": [
            {
                "name": 'test-sky',
                "emission": [2, 2, 2],
                "emission_texture": "test-sky"
            }
        ]
    }
    sunsky_lights = {
        "instances": [],
        "environments": [
            {
                "name": 'test-sunsky',
                "emission": [1, 1, 1],
                "emission_texture": "test-sunsky"
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
        # with open(f'tests/{name}.json', 'wt') as f: json.dump(scene, f, indent=4)
        def write_yaml_objects(f, name):
            if name not in scene: return
            if not scene[name]: return
            f.write(name + ":\n")
            for obj in scene[name]:
                f.write('  - name: ' + obj['name'] + '\n')
                for key, value in obj.items():
                    if key == 'name': continue
                    f.write('    ' + key + ': ' + str(value) + '\n')
        with open(f'tests/{name}.yaml', 'wt') as f:
            write_yaml_objects(f, 'cameras')
            write_yaml_objects(f, 'textures')
            write_yaml_objects(f, 'voltextures')
            write_yaml_objects(f, 'materials')
            write_yaml_objects(f, 'shapes')
            write_yaml_objects(f, 'instances')
            write_yaml_objects(f, 'environments')
    make_test('features1', ['test-bunny', 'test-sphere', 'test-bunny', 'test-sphere', 'test-bunny'], ['test-uvgrid', 'test-plastic-sharp', 'test-metal-rough', 'test-plastic-rough', 'test-metal-sharp'], mixed_lights)
    make_test('materials1', ['test-sphere'], ['test-plastic-sharp', 'test-plastic-rough', 'test-matte', 'test-metal-sharp', 'test-metal-rough'], mixed_lights)
    make_test('materials2', ['test-sphere'], ['test-glass-sharp', 'test-glass-rough', 'test-transparent', 'test-thinglass-sharp', 'test-thinglass-rough'], mixed_lights)
    make_test('materials3', ['test-sphere', 'test-sphere', 'test-sphere-displaced', 'test-sphere', 'test-sphere'], ['test-plastic-sharp-bumped', 'test-plastic-sharp-bumped', 'test-matte-displaced', 'test-metal-sharp-bumped', 'test-metal-sharp-bumped'], mixed_lights)
    make_test('shapes1', ['test-sphere', "test-uvsphere-flipcap", "test-disk", "test-uvcylinder", "test-cube"], ['test-uvgrid'], mixed_lights)
    make_test('shapes2', ['test-cube-subdiv', "test-suzanne-subdiv", "test-teapot", "test-bunny", "test-cube-subdiv"], ['test-uvgrid', 'test-plastic-sharp'], mixed_lights)
    make_test('shapes3', ['test-sphere', "test-hairball1", "test-hairball2", "test-hairball3", "test-sphere", "", "test-hairball-interior", "test-hairball-interior", "test-hairball-interior", ""], ['test-matte', 'test-hair', 'test-hair', 'test-hair', 'test-matte'], mixed_lights, xscales=[ 0.5, 1, 1, 1, 0.5 ])
    make_test('arealights1', ['test-bunny', 'test-sphere', 'test-bunny', 'test-sphere', 'test-bunny'], ['test-uvgrid', 'test-plastic-sharp', 'test-metal-rough', 'test-plastic-rough', 'test-metal-sharp'], area_lights)
    make_test('environments1', ['test-bunny', 'test-sphere', 'test-bunny', 'test-sphere', 'test-bunny'], ['test-uvgrid', 'test-plastic-sharp', 'test-metal-rough', 'test-plastic-rough', 'test-metal-sharp'], sunsky_lights)

cli()