#! /usr/bin/env python3 -B

# build utility for easy development
# complete and unreliable hack used for making it easier to develop
# consider this not to be part of yocto

import click, os, platform, markdown, glob, textwrap

def build(target, dirname, buildtype, cmakeopts=''):
    os.system('mkdir -p build/{dirname}; cd build/{dirname}; cmake ../../ -GNinja -DYOCTO_EMBREE=ON -DCMAKE_BUILD_TYPE={buildtype} {cmakeopts}; cmake --build . {target}'.format(target=target, dirname=dirname, buildtype=buildtype, cmakeopts=cmakeopts))

@click.group()
def run():
    pass

@run.command()
@click.argument('target', required=False, default='')
def release(target=''):
    build(target, 'release', 'Release')

@run.command()
@click.argument('target', required=False, default='')
def debug(target=''):
    build(target, 'debug', 'Debug')

@run.command()
@click.argument('target', required=False, default='')
def nogl(target=''):
    build(target, 'nogl', 'Release', '-DYOCTO_OPENGL=OFF')

@run.command()
def xcode():
    os.system('mkdir -p build/xcode; cd build/xcode; cmake -GXcode -DYOCTO_EXPERIMENTAL=ON -DYOCTO_TOOLS=ON ../../; open yocto-gl.xcodeproj')

@run.command()
def clean():
    os.system('rm -rf bin; rm -rf build')

@run.command()
def format():
    for glob in ['yocto/y*.h', 'yocto/y*.cpp', 'apps/y*.cpp', 'apps/y*.h']:
        os.system('clang-format -i -style=file ' + glob)

@run.command()
def formattests():
    import json
    from collections import OrderedDict
    for filename in glob.glob('tests/*.json'):
        def fix_proc(js):
            return js
        with open(filename, "rt") as f: js = json.load(f, object_pairs_hook=OrderedDict)
        js = fix_proc(js)
        with open(filename, "wt") as f: json.dump(js, f, indent=4)

@run.command()
def tests():
    for ext in ['obj', 'gltf', 'json', 'ybin', 'pbrt']:
        os.system(f'rm -rf tests/{ext}; mkdir tests/{ext}')
    for filename in glob.glob('tests/*.json'):
        print(filename)
        basename = os.path.basename(filename).replace('.json','')
        for ext in ['obj', 'gltf', 'json', 'ybin', 'pbrt']:
            opts = '' if ext != 'ybin' else '--build-bvh'
            os.system(f'mkdir tests/{ext}/{basename}')
            os.system(f'mkdir tests/{ext}/{basename}/textures')
            os.system(f'mkdir tests/{ext}/{basename}/meshes')
            os.system(f'./bin/yscnproc {opts} -o tests/{ext}/{basename}/{basename}.{ext} {filename}')

@run.group()
def scenes():
    pass

def run_cmd(cmd):
    print(f'>>> {cmd}')
    os.system(cmd)

@scenes.command('convert')
@click.option('--format', '-F', default='json', type=click.Choice(['json','obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_convert(dirname, scene='*', format='json'):
    for filename in sorted(glob.glob(f'{dirname}/original/{scene}/*.pbrt')):
        srcdir = os.path.dirname(filename)
        outdir = srcdir.replace('original', f'yocto-{format}')
        if os.path.exists(outdir):
            run_cmd(f'rm -rf {outdir} && mkdir -p {outdir}')
    for filename in sorted(glob.glob(f'{dirname}/original/{scene}/*.pbrt')):
        srcdir = os.path.dirname(filename)
        basename = os.path.basename(filename)
        outdir = srcdir.replace('original', f'yocto-{format}')
        outname = basename.replace('.pbrt',f'.{format}')
        # options = '' if 'ecosys' not in srcdir else '--flipyz'
        options = '--uniform-txt'
        run_cmd(f'mkdir -p {outdir}/models && mkdir -p {outdir}/textures')
        run_cmd(f'../yocto-gl/bin/yscnproc {options} {srcdir}/{basename} -o {outdir}/{outname}')
        if os.path.exists(f'{srcdir}/LICENSE.txt'):
            run_cmd(f'cp -r {srcdir}/LICENSE.txt {outdir}/')

@scenes.command('view')
@click.option('--format', '-F', default='json', type=click.Choice(['json','obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scene_view(dirname, scene='*', format='json'):
    for filename in sorted(glob.glob(f'{dirname}/yocto-{format}/{scene}/*.{format}')):
        run_cmd(f'../yocto-gl/bin/yview --eyelight -D {filename}')

@scenes.command('trace')
@click.option('--format', '-F', default='json', type=click.Choice(['json','obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_trace(dirname, scene='*', format='json'):
    for filename in sorted(glob.glob(f'{dirname}/yocto-{format}/{scene}/*.{format}')):
        run_cmd(f'../yocto-gl/bin/yitrace -D {filename}')

@scenes.command('view_gltf')
@click.option('--mode', '-M', required=False, default='glTF')
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_view_gltf(dirname, scene='*', mode='glTF'):
    format = 'gltf'
    for filename in get_filenames(f'{dirname}/yocto', scene, f'{mode}/*.{format}'):
        run_cmd(f'../yocto-gl/bin/yview --eyelight {filename}')

@scenes.command('trace_gltf')
@click.option('--mode', '-M', required=False, default='glTF')
@click.option('--env', '-E', is_flag=True, required=False, default=False)
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_trace_gltf(dirname, scene='*', mode='glTF', env=False):
    format = 'gltf'
    opts = '--add-skyenv -t environment' if env else '-t eyelight'
    for filename in get_filenames(f'{dirname}/yocto', scene, f'{mode}/*.{format}'):
        run_cmd(f'../yocto-gl/bin/yitrace {opts} {filename}')

@scenes.command('render')
@click.option('--resolution', '-r', default=720, type=int)
@click.option('--samples', '-s', default=64, type=int)
@click.option('--tracer', '-t', default='pathtrace')
@click.option('--progressive', '-p', is_flag=True, default=False, type=bool)
@click.option('--format', '-F', default='obj', type=click.Choice(['obj','gltf']))
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_render(dirname, scene='*', format='obj', samples=64, resolution=720, tracer="pathtrace", progressive=False):
    for filename in sorted(glob.glob(f'{dirname}/yocto-{format}/{scene}/*.{format}')):
        outdir = os.path.dirname(filename)
        run_cmd(f'rm {outdir}/*.hdr && rm {outdir}/*.exr')
        outname = filename.replace(f'.{format}','.hdr')
        save_batch = '--save-batch' if progressive else ''
        run_cmd(f'../yocto-gl/bin/ytrace -r {resolution} -s {samples} -t {tracer} {save_batch} -o {outname} {filename}')

@scenes.command('zip')
@click.argument('dirname', required=True)
@click.argument('scene', required=False, default='*')
def scenes_zip(dirname, scene='*'):
    for filename in glob.glob(f'{dirname}/unedited/{scene}'):
        if not os.path.isdir(filename): continue
        basedir = os.path.dirname(filename)
        basename = os.path.basename(filename)
        print(basedir, dirname)
        run_cmd(f'cd {basedir}; zip -r -X -9 -q {basename}.zip {basename}')

@scenes.command('sync')
def scenes_sync():
    os.system('rsync -avc --delete ./ ../yocto-scenes')    

if __name__ == '__main__':
    run()
