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
    default_scene = {
        "cameras": [
            {
                "uri": "cameras/default.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "uri": "cameras/front.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, 0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, 1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "uri": "cameras/back.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.012,
                "focus_distance": distance([0, 0.575, -0.14], [0, 0.05, 0]),
                "frame": lookat([0, 0.575, -1.4], [0, 0.05, 0], [0,1,0])
            },
            {
                "uri": "cameras/perspective-sharp.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "uri": "cameras/perspective-dof.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.025,
                "film_width": 0.036,
                "film_height": 0.015,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "uri": "cameras/orthographic-sharp.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.0,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": True,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
            {
                "uri": "cameras/orthographic-dof.yaml",
                "focal_length": 0.05,
                "lens_aperture": 0.02,
                "film_width": 0.036,
                "film_height": 0.015,
                "orthographic": True,
                "focus_distance": distance([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05]),
                "frame": lookat([-0.75, 0.4, 0.9], [-0.075, 0.05, -0.05], [0,1,0])
            },
        ],
        "textures": [
            {
                "uri": "::yocto::test-floor::textures/test-floor.png"
            },
            {
                "uri": "::yocto::test-uvgrid::textures/test-uvgrid.png"
            },
            {
                "uri": "::yocto::test-bumps::textures/test-bumps.png"
            },
            {
                "uri": "::yocto::test-bumps-normal::textures/test-bumps-normal.png"
            },
            {
                "uri": "::yocto::test-fbm-displacement::textures/test-fbm-displacement.png"
            },
            {
                "uri": "::yocto::test-sky::textures/test-sky.hdr"
            },
            {
                "uri": "::yocto::test-sunsky::textures/test-sunsky.hdr"
            }
        ],
        "materials": [
            {
                "uri": "materials/test-floor.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_texture": "textures/test-floor.png"
            },
            {
                "uri": "materials/test-uvgrid.yaml",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 1, 1, 1 ],
                "roughness": 0.1,
                "diffuse_texture": "textures/test-uvgrid.png"
            },
            {
                "uri": "materials/test-matte.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 1
            },
            {
                "uri": "materials/test-plastic-sharp.yaml",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "roughness": 0.01
            },
            {
                "uri": "materials/test-plastic-rough.yaml",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "roughness": 0.2
            },
            {
                "uri": "materials/test-metal-sharp.yaml",
                "metallic": 1,
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 0
            },
            {
                "uri": "materials/test-metal-rough.yaml",
                "metallic": 1,
                "diffuse": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.2
            },
            {
                "uri": "materials/test-plastic-sharp-bumped.yaml",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "roughness": 0.01,
                "normal_texture": "textures/test-bumps-normal.png"
            },
            {
                "uri": "materials/test-metal-sharp-bumped.yaml",
                "metallic": 1,
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 0,
                "normal_texture": "textures/test-bumps-normal.png"
            },
            {
                "uri": "materials/test-plastic-rough-coated.yaml",
                "specular": [0.04, 0.04, 0.04],
                "coat": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "roughness": 0.2
            },
            {
                "uri": "materials/test-metal-rough-coated.yaml",
                "metallic": 1,
                "coat": [0.04, 0.04, 0.04],
                "diffuse": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.2
            },
            {
                "uri": "materials/test-transparent.yaml",
                "diffuse": [ 0.7, 0.5, 0.5 ],
                "roughness": 1,
                "opacity": 0.2
            },
            {
                "uri": "materials/test-glass-sharp.yaml",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [1, 1, 1],
                "roughness": 0,
                "thin": False
            },
            {
                "uri": "materials/test-glass-rough.yaml",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.1,
                "thin": False
            },
            {
                "uri": "materials/test-thinglass-sharp.yaml",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [1, 1, 1],
                "roughness": 0,
                "thin": True
            },
            {
                "uri": "materials/test-thinglass-rough.yaml",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.05,
                "thin": True
            },
            {
                "uri": "materials/test-hair.yaml",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 1
            },
            {
                "uri": "materials/test-volume-jade.yaml",
                "specular": [0.04, 0.04, 0.04],
                "roughness": 0,
                "transmission": [1, 1, 1],
                "voltransmission": [0.5, 0.5, 0.5],
                "volscatter": [0.3, 0.6, 0.3]
            },
            {
                "uri": "materials/test-volume-cloud.yaml",
                "transmission": [1, 1, 1],
                "voltransmission": [0.65, 0.65, 0.65],
                "volscatter": [0.9, 0.9, 0.9]
            },
            {
                "uri": "materials/test-volume-glass.yaml",
                "specular": [0.04, 0.04, 0.04],
                "roughness": 0,
                "transmission": [1, 1, 1],
                "voltransmission": [1, 0.5, 0.5],
                "volscale": 0.02
            },
            {
                "uri": "materials/test-volume-smoke.yaml",
                "transmission": [1, 1, 1],
                "voltransmission": [0.5, 0.5, 0.5],
                "volscatter": [0.2, 0.2, 0.2],
                "volanisotropy": -0.8
            },
            {
                "uri": "materials/test-volume-emissive.yaml",
                "transmission": [1, 1, 1],
                "voltransmission": [0.95, 0.95, 0.95],
                "volemission": [15, 15, 10],
                "volscatter": [0.01, 0.01, 0.01]
            },
            {
                "uri": "materials/test-arealight1.yaml",
                "emission": [20, 20, 20]
            },
            {
                "uri": "materials/test-arealight2.yaml",
                "emission": [20, 20, 20]
            }
        ],
        "shapes": [
            {
                "uri": "::yocto::test-floor::shapes/test-floor.ply"
            },
            {
                "uri": "shapes/test-bunny.obj"
            },
            {
                "uri": "shapes/test-teapot.obj"
            },
            {
                "uri": "::yocto::test-sphere::shapes/test-sphere.ply"
            },
            {
                "uri": "::yocto::test-cube::shapes/test-cube.ply"
            },
            {
                "uri": "::yocto::test-disk::shapes/test-disk.ply"
            },
            {
                "uri": "::yocto::test-uvsphere-flipcap::shapes/test-uvsphere-flipcap.ply"
            },
            {
                "uri": "::yocto::test-uvcylinder::shapes/test-uvcylinder.ply"
            },
            {
                "uri": "::yocto::test-sphere-displaced::shapes/test-sphere-displaced.obj"
            },
            {
                "uri": "::yocto::test-cube-subdiv::shapes/test-cube-subdiv.obj"
            },
            {
                "uri": "::yocto::test-suzanne-subdiv::shapes/test-suzanne-subdiv.obj"
            },
            {
                "uri": "::yocto::test-hairball1::shapes/test-hairball1.ply"
            },
            {
                "uri": "::yocto::test-hairball2::shapes/test-hairball2.ply"
            },
            {
                "uri": "::yocto::test-hairball3::shapes/test-hairball3.ply"
            },
            {
                "uri": "::yocto::test-hairball-interior::shapes/test-hairball-interior.ply"
            },
            {
                "uri": "::yocto::test-arealight1::shapes/test-arealight1.ply"
            },
            {
                "uri": "::yocto::test-arealight2::shapes/test-arealight2.ply"
            }
        ],
        "subdivs": [
            {
                "uri": "::yocto::test-sphere-displaced::subdivs/test-sphere-displaced.obj",
                "tesselated_shape": "shapes/test-sphere-displaced.obj",
                "preserve_facevarying": True,
                "displacement_texture": "textures/test-fbm-displacement.png",
                "displacement_scale": 0.025
            },
            {
                "uri": "::yocto::test-cube-subdiv::subdivs/test-cube-subdiv.obj",
                "tesselated_shape": "shapes/test-cube-subdiv.obj",
                "subdivision_level": 4,
                "catmull_clark": True,
                "compute_normals": True,
                "preserve_facevarying": True
            },
            {
                "uri": "::yocto::test-suzanne-subdiv::subdivs/test-suzanne-subdiv.obj",
                "tesselated_shape": "shapes/test-suzanne-subdiv.obj",
                "subdivision_level": 2,
                "catmull_clark": True,
                "compute_normals": True
            },
        ],
        "instances": [
            {
                "uri": "instances/test-floor.yaml",
                "shape": "shapes/test-floor.ply",
                "material": "materials/test-floor.yaml"
            }
        ],
        "environments": []
    }
    area_lights = {
        "instances": [
            {
                "uri": "instances/test-arealight1.yaml",
                "shape": "shapes/test-arealight1.ply",
                "material": "materials/test-arealight1.yaml",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "uri": "instances/test-arealight2.yaml",
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
                "uri": "instances/test-arealight1.yaml",
                "shape": "shapes/test-arealight1.ply",
                "material": "materials/test-arealight1.yaml",
                "frame": lookat([ -0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            },
            {
                "uri": "instances/test-arealight2.yaml",
                "shape": "shapes/test-arealight2.ply",
                "material": "materials/test-arealight2.yaml",
                "frame": lookat([ 0.4, 0.8, 0.8 ], [ 0, 0.1, 0 ], [0, 1, 0], True)
            }
        ],
        "environments": [
            {
                "uri": 'environments/test-sky.yaml',
                "emission": [2, 2, 2],
                "emission_texture": "textures/test-sky.hdr"
            }
        ]
    }
    sunsky_lights = {
        "instances": [],
        "environments": [
            {
                "uri": 'environments/test-sunsky.yaml',
                "emission": [2, 2, 2],
                "emission_texture": "textures/test-sunsky.hdr"
            }
        ]
    }
    def make_test(name, shapes, materials, lights, xoffsets=[ -0.4, -0.2, 0, 0.2, 0.4 ], yoffsets=[0,0,0,0,0], zoffsets=[0,0,0,0,0], xscales=[1,1,1,1,1], yscales=[1,1,1,1,1], zscales=[1,1,1,1,1], subdivs=[]):
        import copy
        def remove_preset(filename):
            splits = filename.rpartition('::')
            return splits[2] if splits[2] else splits[0]
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
                'uri': f'instances/{sref}_{mref}.yaml',
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
                    if instance['material'] == remove_preset(material['uri']): used = True
                if used: scene['materials'] += [material] 
            old_subdivs = scene['subdivs']
            scene['subdivs'] = []
            for subdiv in old_subdivs:
                used = False
                if remove_preset(subdiv['uri']) in subdivs:
                    scene['subdivs'] += [subdiv] 
            old_shapes = scene['shapes']
            scene['shapes'] = []
            for shape in old_shapes:
                used = False
                for instance in scene['instances']:
                    if instance['shape'] == remove_preset(shape['uri']): used = True
                for subdiv in scene['subdivs']:
                    if subdiv['tesselated_shape'] == remove_preset(shape['uri']): used = True
                if used: scene['shapes'] += [shape] 
            old_textures = scene['textures']
            scene['textures'] = []
            for texture in old_textures:
                used = False
                for material in scene['materials']:
                    if 'emission_texture' in material and material['emission_texture'] == remove_preset(texture['uri']): used = True
                    if 'diffuse_texture' in material and material['diffuse_texture'] == remove_preset(texture['uri']): used = True
                    if 'normal_texture' in material and material['normal_texture'] == remove_preset(texture['uri']): used = True
                    if 'displacement_texture' in material and material['displacement_texture'] == remove_preset(texture['uri']): used = True
                for subdiv in scene['subdivs']:
                    if 'displacement_texture' in subdiv and subdiv['displacement_texture'] == remove_preset(texture['uri']): used = True
                for environment in scene['environments']:
                    if environment['emission_texture'] == remove_preset(texture['uri']): used = True
                if used: scene['textures'] += [texture] 
        # with open(f'tests/{name}.json', 'wt') as f: json.dump(scene, f, indent=4)
        def write_yaml_objects(f, name):
            if name not in scene: return
            if not scene[name]: return
            f.write(name + ":\n")
            for obj in scene[name]:
                f.write('  - uri: ' + obj['uri'] + '\n')
                for key, value in obj.items():
                    if key == 'uri': continue
                    f.write('    ' + key + ': ' + str(value).lower() + '\n')
        with open(f'{name}', 'wt') as f:
            write_yaml_objects(f, 'cameras')
            write_yaml_objects(f, 'textures')
            write_yaml_objects(f, 'voltextures')
            write_yaml_objects(f, 'materials')
            write_yaml_objects(f, 'shapes')
            write_yaml_objects(f, 'subdivs')
            write_yaml_objects(f, 'instances')
            write_yaml_objects(f, 'environments')
    make_test('tests/features1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], mixed_lights)
    make_test('tests/materials1.yaml', ['shapes/test-sphere.ply'], ["materials/test-plastic-sharp.yaml", "materials/test-plastic-rough.yaml", "materials/test-matte.yaml", "materials/test-metal-sharp.yaml", "materials/test-metal-rough.yaml"], mixed_lights)
    make_test('tests/materials2.yaml', ['shapes/test-sphere.ply'], ["materials/test-glass-sharp.yaml", "materials/test-glass-rough.yaml", "materials/test-transparent.yaml", "materials/test-thinglass-sharp.yaml", "materials/test-thinglass-rough.yaml"], mixed_lights)
    make_test('tests/materials3.yaml', ['shapes/test-sphere.ply', 'shapes/test-sphere.ply', 'shapes/test-sphere.ply', 'shapes/test-sphere.ply', 'shapes/test-sphere.ply'], ["materials/test-plastic-sharp-bumped.yaml", "materials/test-plastic-rough-coated.yaml", "materials/test-metal-sharp-bumped.yaml", "materials/test-metal-rough-coated.yaml", "materials/test-metal-rough.yaml"], mixed_lights)
    make_test('tests/materials4.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-volume-cloud.yaml", "materials/test-volume-glass.yaml", "materials/test-volume-jade.yaml", "materials/test-volume-emissive.yaml", "materials/test-volume-smoke.yaml"], mixed_lights)
    make_test('tests/shapes1.yaml', ['shapes/test-sphere.ply', "shapes/test-uvsphere-flipcap.ply", "shapes/test-disk.ply", "shapes/test-uvcylinder.ply", "shapes/test-cube.ply"], ["materials/test-uvgrid.yaml"], mixed_lights)
    make_test('tests/shapes2.yaml', ['shapes/test-cube-subdiv.obj', "shapes/test-suzanne-subdiv.obj", 'shapes/test-sphere-displaced.obj', "shapes/test-bunny.obj", "shapes/test-teapot.obj"], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-matte.yaml", "materials/test-uvgrid.yaml", "materials/test-uvgrid.yaml"], mixed_lights, subdivs=['subdivs/test-cube-subdiv.obj', "subdivs/test-suzanne-subdiv.obj", "subdivs/test-sphere-displaced.obj"])
    make_test('tests/shapes3.yaml', ['shapes/test-sphere.ply', "shapes/test-hairball1.ply", "shapes/test-hairball2.ply", "shapes/test-hairball3.ply", "shapes/test-sphere.ply", "", "shapes/test-hairball-interior.ply", "shapes/test-hairball-interior.ply", "shapes/test-hairball-interior.ply", ""], ["materials/test-matte.yaml", "materials/test-hair.yaml", "materials/test-hair.yaml", "materials/test-hair.yaml", "materials/test-matte.yaml"], mixed_lights, xscales=[ 0.5, 1, 1, 1, 0.5 ])
    make_test('tests/arealights1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], area_lights)
    make_test('tests/environments1.yaml', ['shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj', 'shapes/test-sphere.ply', 'shapes/test-bunny.obj'], ["materials/test-uvgrid.yaml", "materials/test-plastic-sharp.yaml", "materials/test-metal-rough.yaml", "materials/test-plastic-rough.yaml", "materials/test-metal-sharp.yaml"], sunsky_lights)

cli()