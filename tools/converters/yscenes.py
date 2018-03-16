#! /usr/bin/env python3 -B

# This is just a hacky utility to convert scenes from various format.

import click, glob, os

@click.group()
def run():
    pass

@run.group()
def convert():
    pass

@run.group()
def view():
    pass

@run.group()
def trace():
    pass

@run.group()
def render():
    pass

@run.command()
def sync():
    os.system('rsync -avc ./ ../yocto-scenes')

@convert.command('tungsten')
@click.argument('scene', required=False, default='*')
def convert_tungsten(scene='*'):
    for filename in sorted(glob.glob(f'bitterli/tungsten/{scene}/scene.json')):
        print(f'-------------------------------- {filename}')
        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)
        meshdir = dirname.replace('tungsten','mitsuba')
        outdir = dirname.replace('tungsten','yocto')
        outname = basename.replace('.json','.obj')
        os.system(f'rm -rf {outdir} && mkdir -p {outdir}')
        print(f'../yocto-gl/bin/ytungsten {dirname}/{basename} -o {outdir}/{outname} -m {meshdir}')
        os.system(f'../yocto-gl/bin/ytungsten {dirname}/{basename} -o {outdir}/{outname} -m {meshdir}')
        if os.path.exists(f'{dirname}/textures'):
            os.system(f'cp -r {dirname}/textures {outdir}/')

@view.command('tungsten')
@click.argument('scene', required=False, default='*')
def view_tungsten(scene='*'):
    for filename in sorted(glob.glob(f'bitterli/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/yview -D -c {filename}'
        print(cmd)
        os.system(cmd)

@trace.command('tungsten')
@click.argument('scene', required=False, default='*')
def trace_tungsten(scene='*'):
    for filename in sorted(glob.glob(f'bitterli/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/yitrace -D {filename}'
        print(cmd)
        os.system(cmd)

@render.command('tungsten')
@click.option('--resolution', '-r', default=256, type=int)
@click.option('--samples', '-s', default=64, type=int)
@click.argument('scene', required=False, default='*')
def render_tungsten(scene='*',samples=16,resolution=256):
    for filename in sorted(glob.glob(f'bitterli/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/ytrace -D -r {resolution} -s {samples} {filename}'
        print(cmd)
        os.system(cmd)

@view.command('mcguire')
@click.argument('scene', required=False, default='*')
def view_mcguire(scene='*'):
    for filename in sorted(glob.glob(f'mcguire/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/yview -D -c {filename}'
        print(cmd)
        os.system(cmd)

@trace.command('mcguire')
@click.argument('scene', required=False, default='*')
def trace_mcguire(scene='*'):
    for filename in sorted(glob.glob(f'mcguire/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/yitrace -D {filename}'
        print(cmd)
        os.system(cmd)

@render.command('mcguire')
@click.option('--resolution', '-r', default=256, type=int)
@click.option('--samples', '-s', default=64, type=int)
@click.argument('scene', required=False, default='*')
def render_mcguire(scene='*',samples=16,resolution=256):
    for filename in sorted(glob.glob(f'mcguire/yocto/{scene}/*.obj')):
        print(f'-------------------------------- {filename}')
        cmd = f'../yocto-gl/bin/ytrace -D -r {resolution} -s {samples} {filename}'
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    run()
