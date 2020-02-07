#! /usr/bin/env python3 -B

import click, glob, os, sys, math, json

@click.group()
def cli():
    pass

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='yaml')
@click.option('--mode','-m', default='path')
def itrace(directory='mcguire',scene='*',format='yaml',mode='path'):
    modes = {
        'path': '--bvh highquality',
        'embree': '--bvh embree-highquality',
        'embree-compact': '--bvh embree-compact',
        'eyelight': '-t eyelight --bvh highquality',
        'eyelight-quick': '--all-cameras -s 16 -r 1280 -t eyelight --bvh default'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            cmd = f'../yocto-gl/bin/yscnitrace {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='yaml')
@click.option('--mode','-m', default='default')
def view(directory='mcguire',scene='*',format='yaml',mode='path'):
    modes = {
        'default': '--double-sided',
        'double-sided': '--double-sided',
        'eyelight': '--double-sided --eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            cmd = f'../yocto-gl/bin/yscnview {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='yaml')
@click.option('--mode','-m', default='path')
def trace(directory='mcguire',scene='*',format='yaml',mode='path'):
    modes = {
        'path': '-s 64 -r 640 --bvh highquality',
        'embree': '-s 256 -r 1280 --bvh embree-highquality',
        'embree-compact': '-s 256 -r 1280 --bvh embree-compact',
        'eyelight': '-s 16 -r 1280 -t eyelight --bvh-high-quality',
        'embree-face': '-s 256 -r 640 --bvh embree-highquality',
        'final': '-s 4096 -r 1280 --bvh embree-highquality',
        'final-compact': '-s 4096 -r 1280 --bvh embree-compact',
        'final-filter': '-s 4096 -r 1280 --filter --bvh embree-highquality',
        'final-face': '-s 4096 -r 640 --bvh embree-highquality',
    }
    options = modes[mode]
    outformat = 'png' if 'eyelight' in mode else 'hdr'
    outprefix = 'eyelight' if 'eyelight' in mode else 'images'
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        extracams = []
        if 'sanmiguel' in dirname: extracams = [1, 2]
        if 'island' in dirname: extracams = [1, 2, 3, 4, 5, 6]
        if 'landscape' in dirname: extracams = [1, 2, 3]
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            basename = os.path.basename(filename).replace(f'.{format}','')
            os.system(f'mkdir -p {directory}/{outprefix}-{format}')
            imagename = f'{directory}/{outprefix}-{format}/{basename}.{outformat}'
            cmd = f'../yocto-gl/bin/yscntrace -o {imagename} {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)
            for cam in extracams:
                imagename = f'{directory}/{outprefix}-{format}/{basename}-c{cam}.{outformat}'
                cmd = f'../yocto-gl/bin/yscntrace -o {imagename} --camera {cam} {options} {filename}'
                print(cmd, file=sys.stderr)
                os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='yaml')
@click.option('--mode','-m', default='linear')
def tonemap(directory='mcguire',scene='*',format='yaml',mode='filmic'):
    modes = {
        'linear': '-t --logo',
        'contrast1': '-t --logcontrast 0.6 --logo',
    }
    options = modes[mode]
    outformat = 'png'
    outprefix = 'images'
    from PIL import Image
    from PIL import ImageFont
    from PIL import ImageDraw 
    font = ImageFont.truetype('~/Library/Fonts/FiraSansCondensed-Regular.otf', 18)
    for filename in sorted(glob.glob(f'{directory}/{outprefix}-{format}/{scene}.hdr')+
                           glob.glob(f'{directory}/{outprefix}-{format}/{scene}.exr')):
        imagename = filename.replace(f'.exr',f'.{outformat}').replace(f'.hdr',f'.{outformat}')
        cmd = f'../yocto-gl/bin/yimgproc -o {imagename} {options} {filename}'
        print(cmd, file=sys.stderr)
        os.system(cmd)
        if directory not in ['bitterli', 'disney', 'mcguire', 'pbrt']: continue
        authorfilename = filename.replace('images-yaml/','source/').replace('-fr.','.').replace('-hr.','.').replace('-c1.','.').replace('-c2.','.').replace('-c3.','.').replace('-c4.','.').replace('-c5.','.').replace('-c6.','.').replace('.hdr','') + '/AUTHOR.txt'
        print(authorfilename)
        with open(authorfilename) as f: text = f.read().strip()
        img = Image.open(imagename)
        _, h = img.size
        draw = ImageDraw.Draw(img)
        tw, _ = draw.textsize(text, font=font)
        draw.rectangle([8,h-26-8,8+8+tw,h-8], (0,0,0))
        draw.text((8+4, h-20-8-4),text,(255,255,255),font=font)
        img.save(imagename)
            

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
# @click.option('--mode','-m', default='no-clear')
@click.option('--clean/--no-clean','-C', default=False)
def sync_images(directory='mcguire',scene='*',format='obj',mode='path',clean=True):
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            basename = os.path.basename(filename).replace(f'.{format}','')
            os.system(f'mkdir -p {directory}/images-{format}')
            imagename = f'{directory}/images-{format}/ytrace-{mode}-{basename}.*'
            if clean:
                cmd = f'rm {dirname}/*.png'
                print(cmd, file=sys.stderr)
                os.system(cmd)
                cmd = f'rm {dirname}/*.hdr'
                print(cmd, file=sys.stderr)
                os.system(cmd)
            cmd = f'cp {imagename} {dirname}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--outformat','-F', default='yaml')
@click.option('--mode','-m', default='default')
@click.option('--clean/--no-clean','-C', default=False)
def convert(directory='mcguire',scene='*',format='obj',outformat="yaml",mode='path',clean=True):
    modes = {
        # 'default': '--uniform-textures --mesh-filenames',
        # 'gltf': '--uniform-textures --mesh-filenames --mesh-directory gltf_meshes/'
        'default': '',
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/source/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        obj_options = ''
        if 'bunny2' in dirname and outformat == 'obj': obj_options = '--obj-instances'
        if 'ecosys' in dirname and outformat == 'obj': obj_options = '--obj-instances'
        if 'landscape' in dirname and outformat == 'obj': obj_options = '--obj-instances'
        if 'fractal' in dirname and outformat == 'obj': obj_options = '--obj-instances'
        if 'pavilion' in dirname and outformat == 'obj': continue
        if 'sanmiguel' in dirname and 'pbrt' in dirname and outformat == 'obj': continue
        outdirname = dirname.replace(f'/source/',f'/{outformat}/')
        if clean: os.system(f'rm -rf {outdirname}')
        os.system(f'mkdir -p {outdirname}')
        os.system(f'mkdir -p {outdirname}/textures')
        if outformat == 'yaml':        
            os.system(f'mkdir -p {outdirname}/shapes')
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            outname = filename.replace(f'/source/',f'/{outformat}/').replace(f'.{format}',f'.{outformat}')
            if format != 'dijson':
                cmd = f'../yocto-gl/bin/yscnproc -o {outname} {options} {obj_options} {filename}'
                print(cmd, file=sys.stderr)
                os.system(cmd)
            else:
                cmd = f'../yocto-gl/bin/yislandproc -o {outname} {options} {obj_options} {filename}'
                print(cmd, file=sys.stderr)
                os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='yuksel')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='hair')
@click.option('--outformat','-F', default='ply')
@click.option('--mode','-m', default='default')
@click.option('--clean-models/--no-clean-models','-C', default=False)
def convert_hair(directory='yuksel',scene='*',format='hair',outformat="ply",mode='path',clean_models=True):
    modes = {
        'default': ''
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        if 'ecosys' in dirname and outformat == 'obj': continue
        if 'landscape' in dirname and outformat == 'obj': continue
        if 'fractal' in dirname and outformat == 'obj': continue
        if 'pavilion' in dirname and outformat == 'obj': continue
        for filename in sorted(glob.glob(f'{dirname}/{format}/*.{format}')):
            outname = filename.replace(f'/{format}/',f'/json/').replace(f'.{format}',f'.{outformat}')
            filedir = os.path.dirname(filename)
            cmd = f'../yocto-gl/bin/ymshproc -o {outname} {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='default')
def backup(directory='mcguire',scene='*',format='obj',mode='default'):
    modes = {
        'default': '-r -X -q',
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        outdir = f'{directory}/backup-{format}'
        basedir = f'{directory}/{format}'
        os.system(f'mkdir -p {outdir}')
        dirname = dirname.replace(basedir+'/','')
        outname = dirname + '.zip'
        os.system(f'rm {outdir}/{outname}')
        cmd = f'cd {basedir}; zip {options} {outname} {dirname}; mv {outname} ../../{outdir}/'
        print(cmd)
        os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='procedurals')
@click.option('--mode','-m', default='skies')
@click.option('--clean/--no-clean','-C', default=False)
def make_procedurals(directory='procedurals',mode='skies',clean=False):
    if mode == 'skies':
        dirname = f'{directory}/hdr/textures'
        os.system(f'mkdir -p {dirname}')
        angles = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90]
        for name in ['sky', 'sun']:
            for angle in angles:
                jsonname = f'{dirname}/_proc.json' 
                outname = f'{dirname}/{name}-{angle:02}.hdr'
                js = {
                    'type': 'sky',
                    'width': 2048,
                    'height': 1024,
                    'sun_angle': math.radians(angle),
                    'has_sun': 'sun' in name,
                    'turbidity': 3
                }
                with open(jsonname, 'w') as f: json.dump(js, f, indent=2)
                cmd = f'../yocto-gl/bin/yimproc -o {outname} {jsonname}'
                print(cmd, file=sys.stderr)
                os.system(cmd)
                os.system(f'rm {jsonname}')
    else:
        print('unknown mode')

@cli.command()
def sync():
    os.system("rsync -avcm --delete --include '*/' --include '*.zip' --include '*.tgz' --include '*.pdf' --exclude='*' ./ ../yocto-scenes")
    # os.system('rsync -avc --delete ./ ../yocto-scenes')

@cli.command()
@click.option('--directory', '-d', default='pbrt-v3-scenes')
@click.option('--scene', '-s', default='*')
def pbrtparse(directory='pbrt-v3-scenes',scene='*'):
    broken_scenes = [
        'bunny-fur/f3-15.pbrt',
        "dambreak/dambreak0.pbrt",
        "hair/curly-hair.pbrt",
        "contemporary-bathroom/contemporary-bathroom.pbrt",
        "head/head.pbrt",
        "ecosys/ecosys.pbrt",
        "sanmiguel/sanmiguel.pbrt",
        "sssdragon/dragon_50.pbrt",
        "white-room/whiteroom-daytime.pbrt",
    ]
    scenes = [
        'barcelona-pavilion/pavilion-day.pbrt',
        'barcelona-pavilion/pavilion-night.pbrt',
        'bathroom/bathroom.pbrt',
        'bmw-m6/bmw-m6.pbrt',
        'breakfast/breakfast.pbrt',
        'buddha-fractal/buddha-fractal.pbrt',
        'bunny-fur/f3-15.pbrt',
        'caustic-glass/glass.pbrt',
        "chopper-titan/chopper-titan.pbrt",
        "cloud/cloud.pbrt",
        "coffee-splash/splash.pbrt",
        "contemporary-bathroom/contemporary-bathroom.pbrt",
        "crown/crown.pbrt",
        "dambreak/dambreak0.pbrt",
        "dragon/f8-4a.pbrt",
        "ecosys/ecosys.pbrt",
        "ganesha/ganesha.pbrt",
        "hair/curly-hair.pbrt",
        "hair/sphere-hairblock.pbr",
        "head/head.pbrt",
        "killeroos/killeroo-simple.pbrt",
        "landscape/view-0.pbrt",
        "lte-orb/lte-orb-silver.pbrt",
        "measure-one/frame85.pbrt",
        "pbrt-book/book.pbrt",
        "sanmiguel/sanmiguel.pbrt",
        "simple/dof-dragons.pbrt",
        "smoke-plume/plume-184.pbrt",
        "sportscar/sportscar.pbrt",
        "sssdragon/dragon_50.pbrt",
        "structuresynth/microcity.pbrt",
        "transparent-machines/frame542.pbrt",
        "tt/tt.pbrt",
        "veach-bidir/bidir.pbrt",
        "veach-mis/mis.pbrt",
        "villa/villa-daylight.pbrt",
        "volume-caustic/caustic.pbrt",
        "vw-van/vw-van.pbrt",
        "white-room/whiteroom-daytime.pbrt",
        "yeahright/yeahright.pbrt",
    ]
    # for filename in scenes:
    #     if scene != '*' and not filename.startswith(f'{scene}/'): continue
    #     cmd = f'../yocto-gl/bin/yitrace {filename}'
    #     print(cmd, file=sys.stderr)
    #     os.system(cmd)
    for filename in scenes:
        if scene != '*' and not filename.startswith(f'{scene}/'): continue
        cmd = f'../yocto-gl/bin/yitrace {directory}/{filename}'
        print(cmd, file=sys.stderr)
        os.system(cmd)

cli()
