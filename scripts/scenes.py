#! /usr/bin/env python3 -B

import click, glob, os, sys, math, json

@click.group()
def cli():
    pass

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='path')
def yitrace(directory='mcguire',scene='*',format='obj',mode='path'):
    modes = {
        'path': '',
        'embree': '--embree',
        'eyelight': '-t eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            cmd = f'../yocto-gl/bin/yitrace {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='default')
def yview(directory='mcguire',scene='*',format='obj',mode='path'):
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
            cmd = f'../yocto-gl/bin/yview {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='path')
def ytrace(directory='mcguire',scene='*',format='obj',mode='path'):
    modes = {
        'path': '-s 64 -r 360',
        'embree': '-s 256 -r 720 --embree',
        'eyelight': '-s 16 -r 720 -t eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        os.system(f'mkdir -p {directory}/{format}/_images')
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            imagename = filename.replace(f'.{format}',f'.{mode}.png')
            imagename = imagename.replace(f'{dirname}',f'{directory}/{format}/_images')
            cmd = f'../yocto-gl/bin/ytrace -o {imagename} {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--outformat','-F', default='json')
@click.option('--mode','-m', default='default')
@click.option('--clean/--no-clean','-C', default=False)
def convert(directory='mcguire',scene='*',format='obj',outformat="json",mode='path',clean=True):
    modes = {
        'default': '--skip-textures --mesh-filenames',
        'gltf': '--skip-textures --mesh-filenames --mesh-directory gltf_meshes/'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{format}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        if '-instanced' in dirname and outformat == 'obj': continue
        outdirname = dirname.replace(f'/{format}/',f'/{outformat}/')
        if clean: os.system(f'rm -rf {outdirname}')
        os.system(f'mkdir -p {outdirname}')
        if outformat != 'obj': os.system(f'mkdir -p {outdirname}/models')
        os.system(f'mkdir -p {outdirname}/textures')
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            if format == 'pbrt':
                with open(filename) as f:
                    if 'WorldBegin' not in f.read(): continue
            outname = filename.replace(f'/{format}/',f'/{outformat}/').replace(f'.{format}',f'.{outformat}')
            cmd = f'../yocto-gl/bin/yscnproc -o {outname} {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)
            cmd = f'cp -r {dirname}/textures {outdirname}'
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
        if '-instanced' in dirname and outformat == 'obj': continue
        for filename in sorted(glob.glob(f'{dirname}/models/*.{format}')):
            outname = filename.replace(f'.{format}',f'.{outformat}')
            filedir = os.path.dirname(filename)
            cmd = f'../yocto-gl/bin/ymshproc -o {outname} {options} {filename}'
            print(cmd, file=sys.stderr)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--mode','-m', default='default')
def zip(directory='mcguire',scene='*',mode='default'):
    modes = {
        'default': '-r -X -q',
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        os.system(f'rm -f {dirname}.zip')
        cmd = f'zip {options} {dirname}.zip {dirname}'
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
