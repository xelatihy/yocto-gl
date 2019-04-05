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
                "name": "cameras/default.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "cameras/front.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, 0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, 1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "name": "cameras/back.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, -0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, -1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "name": "cameras/perspective-sharp.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "cameras/perspective-dof.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.025,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "cameras/orthographic-sharp.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": true,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "name": "cameras/orthographic-dof.yaml",
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
                "name": "textures/test-floor.png.ypreset"
            },
            {
                "name": "textures/test-uvgrid.png.ypreset"
            },
            {
                "name": "textures/test-bump.png.ypreset"
            },
            {
                "name": "textures/test-bump-normal.png.ypreset"
            },
            {
                "name": "textures/test-fbm-displacement.png.ypreset"
            },
            {
                "name": "textures/test-sky.hdr.ypreset"
            },
            {
                "name": "textures/test-sunsky.hdr.ypreset"
            }
        ],
        "materials": [
            {
                "name": "materials/test-floor.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_texture": "textures/test-floor.png"
            },
            {
                "name": "materials/test-uvgrid.yaml",
                "diffuse": [ 1, 1, 1 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1,
                "diffuse_texture": "textures/test-uvgrid.png"
            },
            {
                "name": "materials/test-matte.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
            },
            {
                "name": "materials/test-plastic-sharp.yaml",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01
            },
            {
                "name": "materials/test-plastic-rough.yaml",
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.1
            },
            {
                "name": "materials/test-metal-sharp.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0
            },
            {
                "name": "materials/test-metal-rough.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.1
            },
            {
                "name": "materials/test-matte-displaced.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "displacement_texture": "textures/test-fbm-displacement.png",
                "displacement_scale": 0.025
            },
            {
                "name": "materials/test-plastic-sharp-bumped.yaml",
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "roughness": 0.01,
                "normal_texture": "textures/test-bump-normal.png"
            },
            {
                "name": "materials/test-metal-sharp-bumped.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0,
                "normal_texture": "textures/test-bump-normal.png"
            },
            {
                "name": "materials/test-transparent.yaml",
                "diffuse": [ 0.7, 0.5, 0.5 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1,
                "opacity": 0.2
            },
            {
                "name": "materials/test-glass-sharp.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": true
            },
            {
                "name": "materials/test-glass-rough.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.1,
                "refract": true
            },
            {
                "name": "materials/test-thinglass-sharp.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 1, 1 ],
                "roughness": 0,
                "refract": false
            },
            {
                "name": "materials/test-thinglass-rough.yaml",
                "diffuse": [ 0, 0, 0 ],
                "specular": [ 0.04, 0.04, 0.04 ],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.05,
                "refract": false
            },
            {
                "name": "materials/test-hair.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "specular": [ 0, 0, 0 ],
                "roughness": 1
            },
            {
                "name": "materials/test-arealight1.yaml",
                "emission": [ 20, 20, 20 ]
            },
            {
                "name": "materials/test-arealight2.yaml",
                "emission": [ 20, 20, 20 ]
            }
        ],
        "shapes": [
            {
                "name": "shapes/test-floor.ply.ypreset"
            },
            {
                "name": "shapes/test-bunny.obj"
            },
            {
                "name": "shapes/test-teapot.obj"
            },
            {
                "name": "shapes/test-sphere.ply.ypreset"
            },
            {
                "name": "shapes/test-cube.ply.ypreset"
            },
            {
                "name": "shapes/test-disk.ply.ypreset"
            },
            {
                "name": "shapes/test-uvsphere-flipcap.ply.ypreset"
            },
            {
                "name": "shapes/test-uvcylinder.ply.ypreset"
            },
            {
                "name": "shapes/test-sphere-displaced.obj.ypreset",
                "preserve_facevarying": false
            },
            {
                "name": "shapes/test-cube-subdiv.obj.ypreset",
                "subdivision_level": 4,
                "catmull_clark": true,
                "compute_normals": true,
                "preserve_facevarying": true
            },
            {
                "name": "shapes/test-suzanne-subdiv.obj.ypreset",
                "subdivision_level": 2,
                "catmull_clark": true,
                "compute_normals": true
            },
            {
                "name": "shapes/test-hairball1.ply.ypreset"
            },
            {
                "name": "shapes/test-hairball2.ply.ypreset"
            },
            {
                "name": "shapes/test-hairball3.ply.ypreset"
            },
            {
                "name": "shapes/test-hairball-interior.ply.ypreset"
            },
            {
                "name": "shapes/test-arealight1.ply.ypreset"
            },
            {
                "name": "shapes/test-arealight2.ply.ypreset"
            }
        ],
        "instances": [
            {
                "name": "instances/test-floor.yaml",
                "shape": "shapes/test-floor.ply",
                "material": "materials/test-floor.yaml"
            }
        ],
        "environments": []
    }
    area_lights = {
        "instances": [
            {
                "name": "instances/test-arealight1.yaml",
                "shape": "shapes/test-arealight1.ply",
                "material": "materials/test-arealight1.yaml",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "name": "instances/test-arealight2.yaml",
                "shape": "shapes/test-arealight2.ply",
                "material": "materials/test-arealight2.yaml",
                "frame": lookat([ 0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            }
        ],
        "environments": []
    }
    mixed_lights = {
        "instances": [
            {
                "name": "instances/test-arealight1.yaml",
                "shape": "shapes/test-arealight1.ply",
                "material": "materials/test-arealight1.yaml",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "name": "instances/test-arealight2.yaml",
                "shape": "shapes/test-arealight2.ply",
                "material": "materials/test-arealight2.yaml",
                "frame": lookat([ 0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            }
        ],
        "environments": [
            {
                "name": 'environments/test-sky.yaml',
                "emission": [2, 2, 2],
                "emission_texture": "textures/test-sky.hdr"
            }
        ]
    }
    sunsky_lights = {
        "instances": [],
        "environments": [
            {
                "name": 'environments/test-sunsky.yaml',
                "emission": [1, 1, 1],
                "emission_texture": "textures/test-sunsky.hdr"
            }
        ]
    }
    def make_test(name, shapes, materials, lights, xoffsets=[ -0.4, -0.2, 0, 0.2, 0.4 ], yoffsets=[0,0,0,0,0], zoffsets=[0,0,0,0,0], xscales=[1,1,1,1,1], yscales=[1,1,1,1,1], zscales=[1,1,1,1,1]):
        import copy
        def remove_preset(filename):
            return filename.replace('.ypreset', '')
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
            sref = os.path.splitext(os.path.basename(remove_preset(shape)))[0]
            mref = os.path.splitext(os.path.basename(remove_preset(material)))[0]
            scene['instances'] += [ {
                'name': f'instances/{sref}_{mref}.yaml',
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
                    if instance['material'] == remove_preset(material['name']): used = True
                if used: scene['materials'] += [material] 
            old_shapes = scene['shapes']
            scene['shapes'] = []
            for shape in old_shapes:
                used = False
                for instance in scene['instances']:
                    if instance['shape'] == remove_preset(shape['name']): used = True
                if used: scene['shapes'] += [shape] 
            old_textures = scene['textures']
            scene['textures'] = []
            for texture in old_textures:
                used = False
                for material in scene['materials']:
                    if 'emission_texture' in material and material['emission_texture'] == remove_preset(texture['name']): used = True
                    if 'diffuse_texture' in material and material['diffuse_texture'] == remove_preset(texture['name']): used = True
                    if 'normal_texture' in material and material['normal_texture'] == remove_preset(texture['name']): used = True
                    if 'displacement_texture' in material and material['displacement_texture'] == remove_preset(texture['name']): used = True
                for environment in scene['environments']:
                    if environment['emission_texture'] == remove_preset(texture['name']): used = True
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
        with open(f'{name}', 'wt') as f:
            write_yaml_objects(f, 'cameras')
            write_yaml_objects(f, 'textures')
            write_yaml_objects(f, 'voltextures')
            write_yaml_objects(f, 'materials')
            write_yaml_objects(f, 'shapes')
            write_yaml_objects(f, 'instances')
            write_yaml_objects(f, 'environments')
    make_test('tests/features1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], mixed_lights)
    make_test('tests/materials1.yaml', ['shapes/test-sphere.ply'], ["materials/test-plastic-sharp.yaml", "materials/test-plastic-rough.yaml", "materials/test-matte.yaml", "materials/test-metal-sharp.yaml", "materials/test-metal-rough.yaml"], mixed_lights)
    make_test('tests/materials2.yaml', ['shapes/test-sphere.ply'], ["materials/test-glass-sharp.yaml", "materials/test-glass-rough.yaml", "materials/test-transparent.yaml", "materials/test-thinglass-sharp.yaml", "materials/test-thinglass-rough.yaml"], mixed_lights)
    make_test('tests/materials3.yaml', ['shapes/test-sphere.ply', 'shapes/test-sphere.ply', 'shapes/test-sphere-displaced.obj', 'shapes/test-sphere.ply', 'shapes/test-sphere.ply'], ["materials/test-plastic-sharp-bumped.yaml", "materials/test-plastic-sharp-bumped.yaml", "materials/test-matte-displaced.yaml", "materials/test-metal-sharp-bumped.yaml", "materials/test-metal-sharp-bumped.yaml"], mixed_lights)
    make_test('tests/shapes1.yaml', ['shapes/test-sphere.ply', "shapes/test-uvsphere-flipcap.ply", "shapes/test-disk.ply", "shapes/test-uvcylinder.ply", "shapes/test-cube.ply"], ["materials/test-uvgrid.yaml"], mixed_lights)
    make_test('tests/shapes2.yaml', ['shapes/test-cube-subdiv.obj', "shapes/test-suzanne-subdiv.obj", "shapes/test-teapot.obj", "shapes/test-bunny.obj", "shapes/test-cube-subdiv.obj"], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml"], mixed_lights)
    make_test('tests/shapes3.yaml', ['shapes/test-sphere.ply', "shapes/test-hairball1.ply", "shapes/test-hairball2.ply", "shapes/test-hairball3.ply", "shapes/test-sphere.ply", "", "shapes/test-hairball-interior.ply", "shapes/test-hairball-interior.ply", "shapes/test-hairball-interior.ply", ""], ["materials/test-matte.yaml", "materials/test-hair.yaml", "materials/test-hair.yaml", "materials/test-hair.yaml", "materials/test-matte.yaml"], mixed_lights, xscales=[ 0.5, 1, 1, 1, 0.5 ])
    make_test('tests/arealights1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], area_lights)
    make_test('tests/environments1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], sunsky_lights)

cli()