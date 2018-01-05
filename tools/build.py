#! /usr/bin/env python3 -B

# build utility for easy development
# complete and unreliable hack used for making it easier to develop

import click, os, platform, markdown, glob, textwrap

def build(target, dirname, buildtype, cmakeopts=''):
    os.system('mkdir -p build/{dirname}; cd build/{dirname}; cmake ../../ -GNinja -DCMAKE_BUILD_TYPE={buildtype} {cmakeopts}; cmake --build . {target}'.format(target=target, dirname=dirname, buildtype=buildtype, cmakeopts=cmakeopts))
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
    os.system('mkdir -p build/xcode; cd build/xcode; cmake -G Xcode ../../; open yocto-gl.xcodeproj')

@run.command()
def clean():
    os.system('rm -rf bin; rm -rf build')

@run.command()
def format():
    for glob in ['yocto/yocto_*.h', 'yocto/yocto_*.cpp', 'apps/y*.cpp']:
        os.system('clang-format -i -style=file ' + glob)

@run.command()
def docs():
    os.system('./tools/cpp2doc.py')

@run.command()
def deploy():
    def run(cmd):
        print(cmd)
        os.system(cmd)
    run('cp CMakeLists.txt ../yocto-gl/CMakeLists.txt')
    run('cp readme.md ../yocto-gl/')
    run('cp .travis.yml ../yocto-gl/')
    run('cp appveyor.yml ../yocto-gl/')
    run('cp yocto/yocto_* ../yocto-gl/yocto/')
    run('cp yocto/ext/* ../yocto-gl/yocto/ext/')
    run('cp yocto/ext/imgui/* ../yocto-gl/yocto/ext/imgui/')
    run('cp apps/* ../yocto-gl/apps/')
    run('cp apps/w32/* ../yocto-gl/apps/w32/')
    run('cp apps/w32/include/GL/* ../yocto-gl/apps/w32/include/GL/')
    run('cp apps/w32/include/GLFW/* ../yocto-gl/apps/w32/include/GLFW/')
    run('cp apps/w32/lib-vc2015/* ../yocto-gl/apps/w32/lib-vc2015/')
    run('cp images/* ../yocto-gl/images/')
    run('cp tools/* ../yocto-gl/tools/')
    run('cp tools/gltf/schema/* ../yocto-gl/tools/gltf/schema/')
    run('cp tools/gltf/types/* ../yocto-gl/tools/gltf/types/')
    run('cp docs/* ../yocto-gl/docs/')
    run('cp docs/images/* ../yocto-gl/docs/images/')

@run.command()
def readme():
    readme = ''
    for filename in ['yocto/yocto_gl.h']:
        with open(filename) as f: cpp = f.read();
        for line in cpp.splitlines(True):
            if not line.startswith('///'): break
            readme += line.replace('/// ','').replace('///','')
    with open('readme_.md', 'w') as f: f.write(readme)

if __name__ == '__main__':
    run()
