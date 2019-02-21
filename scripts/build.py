#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
def release():
    os.system('mkdir -p build && mkdir -p build/release && cd build/release && cmake ../.. -GNinja -DYOCTO_EMBREE=ON && ninja')

@cli.command()
def debug():
    os.system('mkdir -p build && mkdir -p build/debug && cd build/debug && cmake ../.. -GNinja -DYOCTO_EMBREE=ON -DCMAKE_BUILD_TYPE=Debug && ninja')

@cli.command()
def xcode():
    os.system('mkdir -p build && mkdir -p build/xcode && cd build/xcode && cmake ../.. -GXcode -DYOCTO_EMBREE=ON && open yocto-gl.xcodeproj')

@cli.command()
def clean():
    os.system('rm -rf bin && rm -rf build')

@cli.command()
def tidy():
    os.system('/usr/local/opt/llvm/bin/clang-tidy -checks="-* readability-*" yocto/*.cpp -- -std=c++17 -I./yocto -I/usr/local/include -I/opt/local/include -I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1 -I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/10.0.0/include -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include')

@cli.command()
def tests():
    for ext in ['obj', 'gltf', 'json', 'ybin', 'pbrt']:
        os.system(f'rm -rf tests/{ext}; mkdir tests/{ext}')
    for filename in glob.glob('tests/*.json'):
        print(filename)
        basename = os.path.basename(filename).replace('.json','')
        for ext in ['obj', 'gltf', 'json', 'ybin', 'pbrt']:
            os.system(f'mkdir -p tests/{ext}/{basename}')
            os.system(f'mkdir -p tests/{ext}/{basename}/textures')
            os.system(f'mkdir -p tests/{ext}/{basename}/meshes')
            os.system(f'mkdir -p tests/{ext}/{basename}/models')
            os.system(f'./bin/yscnproc -o tests/{ext}/{basename}/{basename}.{ext} {filename}')

@cli.command()
def format():
    os.system('clang-format -i -style=file yocto/y*.h')
    os.system('clang-format -i -style=file yocto/y*.cpp')
    os.system('clang-format -i -style=file apps/y*.h')
    os.system('clang-format -i -style=file apps/y*.cpp')

cli()