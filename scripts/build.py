#! /usr/bin/env python3 -B

import click
import glob
import os


@click.group()
def cli():
    pass


@cli.command()
@click.option('--clear/--no-clear', default=False)
def release(clear=False):
    os.makedirs('build/terminal/Release', exist_ok=True)
    os.chdir('build/terminal/Release')
    os.system('cmake ../../.. -GNinja -DCMAKE_BUILD_TYPE=Release -DYOCTO_EMBREE=ON')
    os.system('cmake --build . --parallel 8' +
              (' --clean-first' if clear else ''))


@cli.command()
@click.option('--clear/--no-clear', default=False)
def debug(clear=False):
    os.makedirs('build/terminal/Debug', exist_ok=True)
    os.chdir('build/terminal/Debug')
    os.system('cmake ../../.. -GNinja -DCMAKE_BUILD_TYPE=Debug -DYOCTO_EMBREE=ON')
    os.system('cmake --build . --parallel 8' +
              (' --clean-first' if clear else ''))


@cli.command()
def xcode():
    os.makedirs('build/xcode', exist_ok=True)
    os.chdir('build/xcode')
    os.system('cmake ../.. -GXcode -DYOCTO_EMBREE=ON')
    os.system('open yocto_gl.xcodeproj')


@cli.command()
def vs():
    os.makedirs('build/vs', exist_ok=True)
    os.chdir('build/vs')
    os.system('cmake ../.. -G  "Visual Studio 15 2017" -DYOCTO_EMBREE=ON')
    os.system('yocto_gl.sln')


@cli.command()
def clean():
    os.system('rm -rf bin && rm -rf build')


@cli.command()
def format():
    filenames = sorted(glob.glob('libs/*/y*.h') +
                       glob.glob('libs/*/y*.cpp') + glob.glob('apps/*/y*.cpp'))
    for filename in filenames:
        if 'yshapedata' in filename or 'yshape_data' in filename or 'yscene_data' in filename:
            continue
        print(f'formatting {filename}')
        os.system(f'clang-format -i -style=file {filename}')


cli()
