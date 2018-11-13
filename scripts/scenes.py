#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='path')
def yitrace(directory='mcguire',scene='*',format='obj',mode='path'):
    modes = {
        'path': '',
        'embree': '--embree',
        'eyelight': '-t eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            cmd = f'../yocto-gl/bin/yitrace {options} {filename}'
            print(cmd)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='default')
def yview(directory='mcguire',scene='*',format='obj',mode='path'):
    modes = {
        'default': '',
        'double-sided': '--double-sided',
        'eyelight': '--double-sided --eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            cmd = f'../yocto-gl/bin/yview {options} {filename}'
            print(cmd)
            os.system(cmd)

@cli.command()
@click.option('--directory', '-d', default='mcguire')
@click.option('--scene', '-s', default='*')
@click.option('--format','-f', default='obj')
@click.option('--mode','-m', default='path')
def ytrace(directory='mcguire',scene='*',format='obj',mode='path'):
    modes = {
        'path': '-s 64 -r 360',
        'embree': '-s 64 -r 360',
        'eyelight': '-s 16 -r 720 -t eyelight'
    }
    options = modes[mode]
    for dirname in sorted(glob.glob(f'{directory}/{scene}')):
        if not os.path.isdir(dirname): continue
        if '/_' in dirname: continue
        for filename in sorted(glob.glob(f'{dirname}/*.{format}')):
            imagename = filename.replace(f'.{format}',f'.png')
            imagename = imagename.replace(f'{dirname}/',f'{dirname}/ytrace-{mode}.')
            imagename = imagename.replace(f'{dirname}',f'{directory}/_images')
            cmd = f'../yocto-gl/bin/ytrace -o {imagename} {options} {filename}'
            print(cmd)
            os.system(cmd)

@cli.command()
def sync():
    os.system('rsync -avc --delete ./ ../yocto-scenes')

cli()
