#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
def release():
    os.system('mkdir -p build && mkdir -p build/release && cd build/release && cmake ../.. -GNinja -DYOCTO_EMBREE=ON')

@cli.command()
def debug():
    os.system('mkdir -p build && mkdir -p build/debug && cd build/debug && cmake ../.. -GNinja -DYOCTO_EMBREE=ON -DCMAKE_BUILD_TYPE=Debug')

@cli.command()
def xcode():
    os.system('mkdir -p build && mkdir -p build/xcode && cd build/xcode && cmake ../.. -GXcode -DYOCTO_EMBREE=ON && open yocto-gl.xcodeproj')

@cli.command()
def clean():
    os.system('rm -rf bin && rm -rf build')

@cli.command()
def format():
    os.system('clang-format -i -style=file yocto/y*.h')
    os.system('clang-format -i -style=file yocto/y*.cpp')
    os.system('clang-format -i -style=file apps/y*.h')
    os.system('clang-format -i -style=file apps/y*.cpp')

cli()