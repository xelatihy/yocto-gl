#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
@click.option('--clear/--no-clear', default=False)
def release(clear=False):
    os.makedirs('build/terminal/Release', exist_ok=True)
    os.chdir('build/terminal/Release')
    os.system('cmake ../../.. -GNinja -DCMAKE_BUILD_TYPE=Release -DYOCTO_EMBREE=ON')
    os.system('cmake --build . --parallel 8' + (' --clean-first' if clear else ''))

@cli.command()
@click.option('--clear/--no-clear', default=False)
def debug(clear=False):
    os.makedirs('build/terminal/Debug', exist_ok=True)
    os.chdir('build/terminal/Debug')
    os.system('cmake ../../.. -GNinja -DCMAKE_BUILD_TYPE=Debug -DYOCTO_EMBREE=ON')
    os.system('cmake --build . --parallel 8' + (' --clean-first' if clear else ''))

@cli.command()
def xcode():
    os.makedirs('build/xcode', exist_ok=True)
    os.chdir('build/xcode')
    os.system('cmake ../.. -GXcode -DYOCTO_EMBREE=ON')
    os.system('open yocto-gl.xcodeproj')

@cli.command()
def clean():
    os.system('rm -rf bin && rm -rf build')

@cli.command()
def format():
    os.system('clang-format -i -style=file libs/*/y*.h')
    os.system('clang-format -i -style=file libs/*/y*.cpp')
    os.system('clang-format -i -style=file apps/*/y*.h')
    os.system('clang-format -i -style=file apps/*/y*.cpp')

cli()