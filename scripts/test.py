#! /usr/bin/env python3 -B

import click, glob, os, json, yaml

@click.group()
def cli():
    pass

@cli.command()
@click.option('--scene', '-s', default='*.yaml')
def render(scene='*.yaml'):
    for filename in sorted(glob.glob(f'tests/{scene}')):
        print(f'rendering {filename}')
        imfilename = filename.replace('.yaml','.hdr')
        os.system(f'./bin/yscntrace {filename} -o {imfilename} -s 1024')

@cli.command()
@click.option('--image', '-i', default='*.hdr')
def tonemap(image='*.yaml'):
    from PIL import Image
    from PIL import ImageFont
    from PIL import ImageDraw 
    font = ImageFont.truetype('~/Library/Fonts/FiraSansCondensed-Regular.otf', 18)
    msg = {
        'features1': 'Example materials: matte, plastic, metal, glass, subsurface, normal mapping',
        'features2': 'Example shapes: procedural shapes, Catmull-Clark subdivision, hairs, displacement mapping',
    }
    for filename in sorted(glob.glob(f'tests/{image}')):
        outname = filename.replace('.hdr','.png')
        os.system(f'./bin/yimgproc {filename} -o {outname} -t --logo')
        text = ''
        for k in msg:
            if k in filename: text = msg[k]
        if not text: continue
        img = Image.open(outname)
        _, h = img.size
        draw = ImageDraw.Draw(img)
        tw, _ = draw.textsize(text, font=font)
        draw.rectangle([8,h-26-8,8+8+tw,h-8], (0,0,0))
        draw.text((8+4, h-20-8-4),text,(255,255,255),font=font)
        img.save(outname)

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
                "name": "default",
                "lens": 0.05,
                "aperture": 0.0,
                "aspect": 2.4,
                "lookat": [-0.75, 0.4, 0.9, -0.075, 0.05, -0.05, 0,1,0]
            },
            {
                "name": "front",
                "lens": 0.05,
                "aspect": 3,
                "aperture": 0.0,
                "lookat": [0, 0.575, 1.4, 0, 0.05, 0, 0,1,0]
            },
            {
                "name": "back",
                "lens": 0.05,
                "aperture": 0.0,
                "aspect": 3,
                "lookat": [0, 0.575, -1.4, 0, 0.05, 0, 0,1,0]
            },
            {
                "name": "perspective-sharp",
                "lens": 0.05,
                "aperture": 0.0,
                "aspect": 2.4,
                "lookat": [-0.75, 0.4, 0.9, -0.075, 0.05, -0.05, 0,1,0]
            },
            {
                "name": "perspective-dof",
                "lens": 0.05,
                "aperture": 0.025,
                "aspect": 2.4,
                "lookat": [-0.75, 0.4, 0.9, -0.075, 0.05, -0.05, 0,1,0]
            },
            {
                "name": "orthographic-sharp",
                "lens": 0.05,
                "aperture": 0.0,
                "aspect": 2.4,
                "orthographic": True,
                "lookat": [-0.75, 0.4, 0.9, -0.075, 0.05, -0.05, 0,1,0]
            },
            {
                "name": "orthographic-dof",
                "lens": 0.05,
                "aperture": 0.02,
                "aspect": 2.4,
                "orthographic": True,
                "lookat": [-0.75, 0.4, 0.9, -0.075, 0.05, -0.05, 0,1,0]
            },
        ],
        "cameras1": [
            {
                "name": "default",
                "lens": 0.1,
                "aperture": 0.0,
                "aspect": 2.4,
                "lookat": [-0.6, 1.5, 2.75, -0.05, 0.15, 0, 0,1,0]
            },
            {
                "name": "front",
                "lens": 0.1,
                "aperture": 0.0,
                "aspect": 3,
                "lookat": [0, 1.75, 3, 0, 0.175, 0, 0,1,0]
            },
            {
                "name": "back",
                "lens": 0.1,
                "aperture": 0.0,
                "aspect": 3,
                "lookat": [0, 1.5, -3.25, 0, -0.05, 0, 0,1,0]
            },
            {
                "name": "perspective-sharp",
                "lens": 0.1,
                "aperture": 0.0,
                "aspect": 2.4,
                "lookat": [-0.6, 1.5, 2.75, -0.05, 0.15, 0, 0,1,0]
            },
            {
                "name": "perspective-dof",
                "lens": 0.1,
                "aperture": 0.05,
                "aspect": 2.4,
                "lookat": [-0.6, 1.5, 2.75, -0.05, 0.15, 0, 0,1,0]
            },
            {
                "name": "orthographic-sharp",
                "lens": 0.03,
                "aperture": 0.0,
                "aspect": 2.4,
                "orthographic": True,
                "lookat": [-0.5, 1, 2, -0.05, 0.15, 0, 0,1,0]
            },
            {
                "name": "orthographic-dof",
                "lens": 0.03,
                "aperture": 0.02,
                "aspect": 2.4,
                "orthographic": True,
                "lookat": [-0.5, 1, 2, -0.05, 0.15, 0, 0,1,0]
            },
        ],
        "textures": [
            {
                "name": "floor", "filename": "textures/test-floor.png", "preset": "test-floor"
            },
            {
                "name": "uvgrid", "filename": "textures/test-uvgrid.png", "preset": "test-uvgrid"
            },
            {
                "name": "bumps", "filename": "textures/test-bumps.png", "preset": "test-bumps"
            },
            {
                "name": "bumps-normal", "filename": "textures/test-bumps-normal.png", "preset": "test-bumps-normal"
            },
            {
                "name": "bumps-displacement", "filename": "textures/test-bumps-displacement.png", "preset": "test-bumps-displacement"
            },
            {
                "name": "fbm-displacement", "filename": "textures/test-fbm-displacement.png", "preset": "test-fbm-displacement"
            },
            {
                "name": "sky", "filename": "textures/test-sky.hdr", "preset": "test-sky"
            },
            {
                "name": "sunsky", "filename": "textures/test-sunsky.hdr", "preset": "test-sunsky"
            }
        ],
        "materials": [
            {
                "name": "floor",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "diffuse_tex": "floor"
            },
            {
                "name": "uvgrid",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 1, 1, 1 ],
                "roughness": 0.1,
                "diffuse_tex": "uvgrid"
            },
            {
                "name": "matte",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 1
            },
            {
                "name": "plastic-sharp",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "roughness": 0.01
            },
            {
                "name": "plastic-rough",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "roughness": 0.2
            },
            {
                "name": "metal-sharp",
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0
            },
            {
                "name": "metal-rough",
                "specular": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.2
            },
            {
                "name": "plastic-sharp-bumped",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.5, 0.7 ],
                "roughness": 0.01,
                "normal_tex": "bumps-normal"
            },
            {
                "name": "plastic-rough-bumped",
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "roughness": 0.2,
                "normal_tex": "bumps-normal"
            },
            {
                "name": "metal-sharp-bumped",
                "specular": [ 0.7, 0.7, 0.7 ],
                "roughness": 0,
                "normal_tex": "bumps-normal"
            },
            {
                "name": "plastic-rough-coated",
                "specular": [0.04, 0.04, 0.04],
                "coat": [0.04, 0.04, 0.04],
                "diffuse": [ 0.5, 0.7, 0.5 ],
                "roughness": 0.2
            },
            {
                "name": "metal-rough-coated",
                "coat": [0.04, 0.04, 0.04],
                "specular": [ 0.66, 0.45, 0.34 ],
                "roughness": 0.2
            },
            {
                "name": "uvgrid-coated",
                "coat": [0.04, 0.04, 0.04],
                "specular": [0.04, 0.04, 0.04],
                "diffuse": [ 1, 1, 1 ],
                "roughness": 0.2,
                "diffuse_tex": "uvgrid"
            },
            {
                "name": "transparent",
                "diffuse": [ 0.7, 0.5, 0.5 ],
                "roughness": 1,
                "opacity": 0.2
            },
            {
                "name": "glass-sharp",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [1, 1, 1],
                "roughness": 0,
                "refract": True
            },
            {
                "name": "glass-rough",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.1,
                "refract": True
            },
            {
                "name": "thinglass-sharp",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [1, 1, 1],
                "roughness": 0
            },
            {
                "name": "thinglass-rough",
                "specular": [0.04, 0.04, 0.04],
                "transmission": [ 1, 0.7, 0.7 ],
                "roughness": 0.05
            },
            {
                "name": "hair",
                "diffuse": [ 0.7, 0.7, 0.7 ],
                "roughness": 1
            },
            {
                "name": "volume-jade",
                "specular": [0.04, 0.04, 0.04],
                "roughness": 0,
                "transmission": [1, 1, 1],
                "voltransmission": [0.5, 0.5, 0.5],
                "volscatter": [0.3, 0.6, 0.3],
                "refract": True            
            },
            {
                "name": "volume-cloud",
                "transmission": [1, 1, 1],
                "voltransmission": [0.65, 0.65, 0.65],
                "volscatter": [0.9, 0.9, 0.9]
            },
            {
                "name": "volume-glass",
                "specular": [0.04, 0.04, 0.04],
                "roughness": 0,
                "transmission": [1, 1, 1],
                "voltransmission": [1, 0.5, 0.5],
                "volscale": 0.02,
                "refract": True
            },
            {
                "name": "volume-smoke",
                "transmission": [1, 1, 1],
                "voltransmission": [0.5, 0.5, 0.5],
                "volscatter": [0.2, 0.2, 0.2],
                "volanisotropy": -0.8
            },
            {
                "name": "volume-emissive",
                "transmission": [1, 1, 1],
                "voltransmission": [0.95, 0.95, 0.95],
                "volemission": [15, 15, 10],
                "volscatter": [0.01, 0.01, 0.01]
            },
            {
                "name": "arealight1",
                "emission": [20, 20, 20]
            },
            {
                "name": "arealight2",
                "emission": [20, 20, 20]
            },
            {
                "name": "largearealight1",
                "emission": [10, 10, 10]
            },
            {
                "name": "largearealight2",
                "emission": [10, 10, 10]
            }
        ],
        "shapes": [
            {
                "name": "floor", "filename": "shapes/test-floor.ply", "preset": "test-floor"
            },
            {
                "name": "bunny", "filename": "shapes/test-bunny.obj"
            },
            {
                "name": "teapot", "filename": "shapes/test-teapot.obj"
            },
            {
                "name": "sphere", "filename": "shapes/test-sphere.ply", "preset": "test-sphere"
            },
            {
                "name": "cube", "filename": "shapes/test-cube.ply", "preset": "test-cube"
            },
            {
                "name": "disk", "filename": "shapes/test-disk.ply", "preset": "test-disk"
            },
            {
                "name": "uvsphere-flipcap", "filename": "shapes/test-uvsphere-flipcap.ply", "preset": "test-uvsphere-flipcap"
            },
            {
                "name": "uvcylinder", "filename": "shapes/test-uvcylinder.ply", "preset": "test-uvcylinder"
            },
            {
                "name": "sphere-displaced", "filename": "shapes/test-sphere-displaced.obj", "preset": "test-sphere-displaced",
                "facevarying": True,
                "displacement": 0.025,
                "displacement_tex": "bumps-displacement"
            },
            {
                "name": "cube-subdiv", "filename": "shapes/test-cube-subdiv.obj", "preset": "test-cube-subdiv",
                "subdivisions": 4,
                "catmullclark": True,
                "smooth": True,
                "facevarying": True
            },
            {
                "name": "suzanne-subdiv", "filename": "shapes/test-suzanne-subdiv.obj", "preset": "test-suzanne-subdiv",
                "subdivisions": 2,
                "catmullclark": True,
                "smooth": True
            },
            {
                "name": "hairball1", "filename": "shapes/test-hairball1.ply", "preset": "test-hairball1"
            },
            {
                "name": "hairball2", "filename": "shapes/test-hairball2.ply", "preset": "test-hairball2"
            },
            {
                "name": "hairball3", "filename": "shapes/test-hairball3.ply", "preset": "test-hairball3"
            },
            {
                "name": "hairball-interior", "filename": "shapes/test-hairball-interior.ply", "preset": "test-hairball-interior"
            },
            {
                "name": "arealight1", "filename": "shapes/test-arealight1.ply", "preset": "test-arealight1"
            },
            {
                "name": "arealight2", "filename": "shapes/test-arealight2.ply", "preset": "test-arealight2"
            },
            {
                "name": "largearealight1", "filename": "shapes/test-largearealight1.ply", "preset": "test-largearealight1"
            },
            {
                "name": "largearealight2", "filename": "shapes/test-largearealight2.ply", "preset": "test-largearealight2"
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
                "lookat": [ -0.4, 0.8, 0.8 ,  0, 0.1, 0 , 0, 1, 0]
            },
            {
                "name": "arealight2",
                "shape": "arealight2",
                "material": "arealight2",
                "lookat": [ 0.4, 0.8, 0.8 , 0, 0.1, 0 , 0, 1, 0]
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
                "lookat": [ -0.4, 0.8, 0.8, 0, 0.1, 0, 0, 1, 0]
            },
            {
                "name": "arealight2",
                "shape": "arealight2",
                "material": "arealight2",
                "lookat": [ 0.4, 0.8, 0.8, 0, 0.1, 0, 0, 1, 0]
            }
        ],
        "environments": [
            {
                "name": 'sky',
                "emission": [0.5, 0.5, 0.5],
                "emission_tex": "sky"
            }
        ]
    }
    mixed_lights1 = {
        "instances": [
            {
                "name": "largearealight1",
                "shape": "largearealight1",
                "material": "largearealight1",
                "lookat": [ -0.8, 1.6, 1.6 , 0, 0.1, 0, 0, 1, 0]
            },
            {
                "name": "largearealight2",
                "shape": "largearealight2",
                "material": "largearealight2",
                "lookat": [ 0.8, 1.6, 1.6 , 0, 0.1, 0, 0, 1, 0]
            }
        ],
        "environments": [
            {
                "name": 'sky',
                "emission": [0.5, 0.5, 0.5],
                "emission_tex": "sky"
            }
        ]
    }
    sunsky_lights = {
        "instances": [],
        "environments": [
            {
                "name": 'sunsky',
                "emission": [1, 1, 1],
                "emission_tex": "sunsky"
            }
        ]
    }
    def make_test(name, shapes, materials, lights, xoffsets=[ -0.4, -0.2, 0, 0.2, 0.4 ], yoffsets=[0,0,0,0,0], zoffsets=[0,0,0,0,0], xscales=[1,1,1,1,1], yscales=[1,1,1,1,1], zscales=[1,1,1,1,1], inrow=True, new_cameras=False):
        import copy
        def remove_preset(filename):
            splits = filename.rpartition('::')
            return splits[2] if splits[2] else splits[0]
        scene = copy.deepcopy(default_scene)
        if new_cameras:
            scene['cameras'] = scene['cameras1']
        del scene['cameras1']
        scene['instances'] += copy.deepcopy(lights['instances'])
        scene['environments'] += copy.deepcopy(lights['environments'])
        num = max(len(xoffsets), len(shapes), len(materials))
        for i in range(num):
            shape = shapes[i % len(shapes)]
            material = materials[i % len(materials)]
            xoffset = xoffsets[i % len(xoffsets)]
            yoffset = yoffsets[i % len(yoffsets)]
            zoffset = zoffsets[i % len(zoffsets)] if inrow else zoffsets[i // len(xoffsets)]
            xscale = xscales[i % len(xscales)]
            yscale = yscales[i % len(yscales)]
            zscale = zscales[i % len(zscales)]
            if not shape: continue
            if not material: continue
            sref = os.path.splitext(os.path.basename(remove_preset(shape)))[0]
            mref = os.path.splitext(os.path.basename(remove_preset(material)))[0]
            scene['instances'] += [ {
                'name': f'{sref}_{mref}',
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
                    if 'emission_tex' in material and material['emission_tex'] == remove_preset(texture['name']): used = True
                    if 'diffuse_tex' in material and material['diffuse_tex'] == remove_preset(texture['name']): used = True
                    if 'normal_tex' in material and material['normal_tex'] == remove_preset(texture['name']): used = True
                    if 'displacement_tex' in material and material['displacement_tex'] == remove_preset(texture['name']): used = True
                for shape in scene['shapes']:
                    if 'displacement_tex' in shape and shape['displacement_tex'] == remove_preset(texture['name']): used = True
                for environment in scene['environments']:
                    if environment['emission_tex'] == remove_preset(texture['name']): used = True
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
                    f.write('    ' + key + ': ' + str(value).lower() + '\n')
        with open(f'{name}', 'wt') as f:
            write_yaml_objects(f, 'cameras')
            write_yaml_objects(f, 'textures')
            write_yaml_objects(f, 'voltextures')
            write_yaml_objects(f, 'materials')
            write_yaml_objects(f, 'shapes')
            write_yaml_objects(f, 'instances')
            write_yaml_objects(f, 'environments')
    make_test('tests/features1.yaml', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ["uvgrid-coated", "volume-glass", "volume-jade", "plastic-rough-bumped", "metal-rough"], mixed_lights)
    make_test('tests/features2.yaml', ['sphere', 'suzanne-subdiv', 'hairball1', 'sphere-displaced', 'cube', '', '', 'hairball-interior', '', ''], ["uvgrid", "plastic-rough", "hair", "plastic-rough", "uvgrid", '', '', 'hair', '', ''], mixed_lights)
    make_test('tests/materials1.yaml', ['sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'bunny', 'bunny', 'bunny', 'bunny', 'bunny'], ["plastic-sharp", "plastic-rough", "matte", "metal-sharp", "metal-rough", "plastic-sharp", "plastic-rough", "matte", "metal-sharp", "metal-rough"], mixed_lights1, zoffsets=[ 0, -0.4 ], inrow=False, new_cameras=True)
    make_test('tests/materials2.yaml', ['sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'bunny', 'bunny', 'bunny', 'bunny', 'bunny'], ["glass-sharp", "glass-rough", "transparent", "thinglass-sharp", "thinglass-rough", "glass-sharp", "glass-rough", "transparent", "thinglass-sharp", "thinglass-rough"], mixed_lights1, zoffsets=[ 0, -0.4 ], inrow=False, new_cameras=True)
    make_test('tests/materials3.yaml', ['sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'bunny', 'bunny', 'bunny', 'bunny', 'bunny'], ["plastic-sharp-bumped", "plastic-rough-coated", "metal-sharp-bumped", "metal-rough-coated", "metal-rough", "plastic-sharp-bumped", "plastic-rough-coated", "metal-sharp-bumped", "metal-rough-coated", "metal-rough"], mixed_lights1, zoffsets=[ 0, -0.4 ], inrow=False, new_cameras=True)
    make_test('tests/materials4.yaml', ['sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'bunny', 'bunny', 'bunny', 'bunny', 'bunny'], ["volume-cloud", "volume-glass", "volume-jade", "volume-emissive", "volume-smoke", "volume-cloud", "volume-glass", "volume-jade", "volume-emissive", "volume-smoke"], mixed_lights1, zoffsets=[ 0, -0.4 ], inrow=False, new_cameras=True)
    make_test('tests/shapes1.yaml', ['sphere', "uvsphere-flipcap", "disk", "uvcylinder", "cube"], ["uvgrid"], mixed_lights)
    make_test('tests/shapes2.yaml', ['cube-subdiv', "suzanne-subdiv", 'sphere-displaced', "bunny", "teapot"], ["uvgrid", "plastic-sharp", "matte", "uvgrid", "uvgrid"], mixed_lights)
    make_test('tests/shapes3.yaml', ['sphere', "hairball1", "hairball2", "hairball3", "sphere", "", "hairball-interior", "hairball-interior", "hairball-interior", ""], ["matte", "hair", "hair", "hair", "matte"], mixed_lights, xscales=[ 0.5, 1, 1, 1, 0.5 ])
    make_test('tests/arealights1.yaml', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ["uvgrid", "plastic-sharp", "metal-rough", "plastic-rough", "metal-sharp"], area_lights)
    make_test('tests/environments1.yaml', ['bunny', 'sphere', 'bunny', 'sphere', 'bunny'], ["uvgrid", "plastic-sharp", "metal-rough", "plastic-rough", "metal-sharp"], sunsky_lights)
    make_test('tests/materials.yaml', ['bunny'], 
        ["plastic-sharp", "plastic-rough", "matte", "metal-sharp", "metal-rough"] +
        ["glass-sharp", "glass-rough", "transparent", "thinglass-sharp", "thinglass-rough"] +
        ["plastic-sharp-bumped", "plastic-rough-coated", "metal-sharp-bumped", "metal-rough-coated", "metal-rough"] +
        ["volume-cloud", "volume-glass", "volume-jade", "volume-emissive", "volume-smoke"],
        mixed_lights1, zoffsets=[ 0.2, 0, -0.2, -0.4 ], inrow=False)

@cli.command()
def upgrade():
    for filename in sorted(glob.glob('tests/*.yaml')):
        continue
        print(filename)
        with open(filename) as f: yaml = f.read()
        nyaml = ''
        cur_type = ''
        cur_name = ''
        materials = {}
        shapes = {}
        for line in yaml.splitlines():
            line = line.rstrip()
            if not line:
                nyaml += line + '\n'
            elif line in ['materials:', 'shapes:', 'instances:']:
                cur_type = line
                if line == 'instances:':
                    nyaml += 'shapes:\n'
            elif not line.startswith('  '):
                cur_type = ''
                nyaml += line + '\n'
            elif line.startswith('  ') and not cur_type:
                nyaml += line + '\n'
            elif line.startswith('  ') and cur_type == 'materials:':
                if line.startswith('  - '):
                    cur_name = line.partition(':')[2].strip()
                    materials[cur_name] = ''
                else:
                    materials[cur_name] += line + '\n'
            elif line.startswith('  ') and cur_type == 'shapes:':
                if line.startswith('  - '):
                    cur_name = line.partition(':')[2].strip()
                    shapes[cur_name] = ''
                else:
                    shapes[cur_name] += line.replace('filename:','shape:') + '\n'
            elif line.startswith('  ') and cur_type == 'instances:':
                if line.startswith('  - '):
                    nyaml += '  - name: ' + line.partition(':')[2].strip().replace('_','-') + '\n'
                elif line.startswith('    material:'):
                    nyaml += materials[line.partition(':')[2].strip()]
                elif line.startswith('    shape:'):
                    nyaml += shapes[line.partition(':')[2].strip()]
                else:
                    nyaml += line + '\n'
        with open(filename, 'wt') as f:
            f.write(nyaml)
    for filename in sorted(glob.glob('tests/*.yaml')):
        print(filename)
        with open(filename) as f: yaml = f.read()
        nyaml = ''
        cur_type = ''
        cur_name = ''
        textures = {}
        for line in yaml.splitlines():
            line = line.rstrip()
            if not line:
                nyaml += line + '\n'
            elif line in ['textures:']:
                cur_type = line
            elif not line.startswith('  '):
                cur_type = ''
                nyaml += line + '\n'
            elif line.startswith('  ') and not cur_type:
                if '_tex:' in line:
                    nyaml += line.partition(':')[0] + ': ' + textures[line.partition(':')[2].strip()] + '\n'
                else:
                    nyaml += line + '\n'
            elif line.startswith('  ') and cur_type == 'textures:':
                if line.startswith('  - '):
                    cur_name = line.partition(':')[2].strip()
                    textures[cur_name] = ''
                else:
                    textures[cur_name] = line.partition(':')[2].strip()
        with open(filename, 'wt') as f:
            f.write(nyaml)

@cli.command()
def fix_dict():
    def noname(d):
        import copy
        d = copy.copy(d)
        del d['name']
        return d
    for filename in sorted(glob.glob('tests/*.json')):
        print(filename)
        with open(filename) as f: scene = json.load(f)
        if 'cameras' in scene: 
            scene['cameras'] = { value['name']: noname(value) for value in scene['cameras'] }
        if 'environments' in scene: 
            scene['environments'] = { value['name']: noname(value) for value in scene['environments'] }
        if 'materials' in scene: 
            scene['materials'] = { value['name']: noname(value) for value in scene['materials'] }
        if 'objects' in scene: 
            scene['objects'] = { value['name']: noname(value) for value in scene['objects'] }
        with open(filename, 'w') as f: json.dump(scene, f, indent=2)

cli()
