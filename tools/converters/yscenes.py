#! /usr/bin/env python3 -B

# This is just a hacky utility to convert scenes from various format.

import click, glob, os

def get_filenames(dirname, sceneglob, fileglob):
    excluded = []
    if os.path.exists(f'{dirname}/excluded.txt'):
        with open(f'{dirname}/excluded.txt') as f:
            excluded = [ f'{dirname}/' + l.strip() for l in f.readlines() ]
    filenames = []
    for filename in sorted(glob.glob(f'{dirname}/{sceneglob}/{fileglob}')):
        is_excluded = False
        for exc in excluded:
            if exc in filename:
                is_excluded = True
        if is_excluded:
            print(f'exclude {filename}')
        else:                    
            filenames += [ filename ]
    return filenames

def run_cmd(cmd):
    print(f'>>> {cmd}')
    os.system(cmd)

@click.group()
def run():
    pass

@run.command()
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def convert(dirname, scene='*'):
    for filename in get_filenames(f'{dirname}/original', scene, '*.pbrt'):
        srcdir = os.path.dirname(filename)
        outdir = srcdir.replace('original','yocto')
        if os.path.exists(outdir):
            run_cmd(f'rm -rf {outdir} && mkdir -p {outdir}')
    for filename in get_filenames(f'{dirname}/original', scene, '*.pbrt'):
        srcdir = os.path.dirname(filename)
        basename = os.path.basename(filename)
        outdir = srcdir.replace('original','yocto')
        objname = basename.replace('.pbrt','.obj')
        gltfname = basename.replace('.pbrt','.gltf')
        # options = '' if 'ecosys' not in srcdir else '--flipyz'
        options = ''
        run_cmd(f'mkdir -p {outdir}/models')
        run_cmd(f'../yocto-gl/bin/ypbrt {options} {srcdir}/{basename} -o {outdir}/{objname}')
        run_cmd(f'../yocto-gl/bin/ypbrt {options} {srcdir}/{basename} -o {outdir}/{gltfname}')
        if os.path.exists(f'{srcdir}/textures'):
            run_cmd(f'cp -r {srcdir}/textures {outdir}/')
        if os.path.exists(f'{srcdir}/LICENSE.txt'):
            run_cmd(f'cp -r {srcdir}/LICENSE.txt {outdir}/')
        for imfilename in glob.glob(f'{srcdir}/textures/*.pfm'):
            out_imfilename = imfilename.replace('.pfm','.hdr')
            run_cmd(f'../yocto-gl/bin/yimproc -o {out_imfilename} {imfilename}')

@run.command()
@click.option('--format', '-F', default='obj', type=click.Choice(['obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def view(dirname, scene='*', format='obj'):
    for filename in get_filenames(f'{dirname}/yocto', scene, f'*.{format}'):
        run_cmd(f'../yocto-gl/bin/yview --eyelight {filename}')

@run.command()
@click.option('--format', '-F', default='obj', type=click.Choice(['obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def trace(dirname, scene='*', format='obj'):
    for filename in get_filenames(f'{dirname}/yocto', scene, f'*.{format}'):
        run_cmd(f'../yocto-gl/bin/yitrace {filename}')

@run.command()
@click.option('--resolution', '-r', default=720, type=int)
@click.option('--samples', '-s', default=64, type=int)
@click.option('--tracer', '-t', default='pathtrace')
@click.option('--progressive', '-p', is_flag=True, default=False, type=bool)
@click.option('--format', '-F', default='obj', type=click.Choice(['obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def render(dirname, scene='*', format='obj', samples=64, resolution=720, tracer="pathtrace", progressive=False):
    for filename in get_filenames(f'{dirname}/yocto', scene, f'*.{format}'):
        outdir = os.path.dirname(filename).replace('/yocto','/render')
        if os.path.exists(f'{outdir}'):
            run_cmd(f'rm {outdir}/*.{format}.exr')
        else:
            run_cmd(f'mkdir -p {outdir}')
        outname = filename.replace('/yocto','/render')+'.exr'
        save_batch = '--save-batch' if progressive else ''
        run_cmd(f'../yocto-gl/bin/ytrace -r {resolution} -s {samples} -t {tracer} {save_batch} -o {outname} {filename}')

@run.command()
def sync():
    os.system('rsync -avc --delete ./ ../yocto-scenes')


if __name__ == '__main__':
    run()
