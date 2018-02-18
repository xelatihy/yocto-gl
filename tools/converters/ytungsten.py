#! /usr/bin/env python3 -B

# Convert tungsten example files from https://benedikt-bitterli.me/resources/ to Yocto/GL.
# This is a hack.

import click, glob, os

@click.group()
def run():
    pass

@run.command()
def convert():
    for filename in sorted(glob.glob('bitterli/tungsten/*/*.json')):
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

@run.command()
def view():
    for filename in sorted(glob.glob('bitterli/yocto/*/*.obj')):
        print(f'-------------------------------- {filename}')
        os.system(f'../yocto-gl/bin/yview {filename}')

if __name__ == '__main__':
    run()
