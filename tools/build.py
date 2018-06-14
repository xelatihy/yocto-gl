#! /usr/bin/env python3 -B

# build utility for easy development
# complete and unreliable hack used for making it easier to develop

import click, os, platform, markdown, glob, textwrap

def build(target, dirname, buildtype, cmakeopts=''):
    os.system('mkdir -p build/{dirname}; cd build/{dirname}; cmake ../../ -GNinja -DCMAKE_BUILD_TYPE={buildtype} -DYOCTO_EXPERIMENTAL=ON -DYOCTO_TOOLS=ON {cmakeopts}; cmake --build . {target}'.format(target=target, dirname=dirname, buildtype=buildtype, cmakeopts=cmakeopts))
    os.system('ln -Ffs {dirname} build/latest'.format(dirname=dirname))

@click.group()
def run():
    pass

@run.command()
@click.argument('target', required=False, default='')
def latest(target=''):
    os.system('cd build/latest; cmake --build . {target}'.format(target=target))

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
@click.argument('target', required=False, default='')
def debug(target=''):
    build(target, 'debug', 'Debug')

@run.command()
@click.argument('target', required=False, default='')
def gcc(target=''):
    build(target, 'gcc', 'Release', '-DCMAKE_C_COMPILER=gcc-7 -DCMAKE_CXX_COMPILER=g++-7')

@run.command()
def xcode():
    os.system('mkdir -p build/xcode; cd build/xcode; cmake -GXcode -DYOCTO_EXPERIMENTAL=ON -DYOCTO_TOOLS=ON ../../; open yocto-gl.xcodeproj')

@run.command()
def clean():
    os.system('rm -rf bin; rm -rf build')

@run.command()
def format():
    for glob in ['yocto/y*.h', 'yocto/y*.cpp', 'apps/y*.cpp']:
        os.system('clang-format -i -style=file ' + glob)

@run.command()
def tests():
    for ext in ['obj', 'gltf', 'json']:
        os.system(f'rm -rf tests/{ext}; mkdir tests/{ext}')
    for filename in glob.glob('tests/*.json'):
        print(filename)
        basename = os.path.basename(filename).replace('.json','')
        for ext in ['obj', 'gltf', 'json']:
            os.system(f'mkdir tests/{ext}/{basename}')
            os.system(f'mkdir tests/{ext}/{basename}/textures')
            os.system(f'mkdir tests/{ext}/{basename}/meshes')
            os.system(f'./bin/yscnproc -o tests/{ext}/{basename}/{basename}.{ext} {filename}')

if __name__ == '__main__':
    run()
