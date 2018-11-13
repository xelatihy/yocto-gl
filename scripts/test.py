#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
def run():
    os.system('mkdir -p build && mkdir -p build/release && cd build/release && cmake ../.. -GNinja -DYOCTO_EMBREE=ON')
    os.system('rm tests/run-tests/output/*.png; rm tests/run-tests/difference/*.png')
    os.system('cd build/release && ctest -j 4 --output-on-failure')

@cli.command()
def clean():
    os.system('rm tests/run-tests/output/*.png; rm tests/run-tests/difference/*.png')

@cli.command()
def update():
    os.system('cp tests/run-tests/output/*.png tests/run-tests/result/')

cli()