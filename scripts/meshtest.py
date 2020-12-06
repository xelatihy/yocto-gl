#! /usr/bin/env python3 -B

import click
import glob
import os
import subprocess
import json

mesh_dir = 'meshes'
scene_dir = 'scenes'
curve_dir = 'curves'
stats_dir = 'stats'


@click.group()
def run():
    pass


@run.command()
@click.option('--dirname', '-d')
def trace(dirname):
    def handle_error(err, result, mesh_name, stats_name):
        msg = 'error: ' + err
        print(msg + ' ' * max(0, 78-len(msg)))
        if err not in result:
            result[err] = []
        result[err] += [mesh_name]
        stats = {'error': err}
        with open(stats_name, 'wt') as f:
            json.dump(stats, f, indent=2)

    result = {'ok': []}
    mesh_names = glob.glob(f'{dirname}/{mesh_dir}/*.obj')
    mesh_num = len(mesh_names)
    for mesh_id, mesh_name in enumerate(mesh_names):
        stats_name = mesh_name.replace(
            'meshes/', 'stats/').replace('.obj', '.json')
        curve_name = mesh_name.replace(
            'meshes/', 'curves/').replace('.obj', '.json')
        scene_name = mesh_name.replace(
            'meshes/', 'scenes/').replace('.obj', '.json')
        msg = f'[{mesh_id}/{mesh_num}] {mesh_name}'
        print(msg + ' ' * max(0, 78-len(msg)))
        cmd = f'../yocto-gl/bin/ymeshtest {mesh_name} -s {stats_name} -S {scene_name} -p {curve_name}'
        try:
            retcode = subprocess.run(cmd, timeout=5, shell=True).returncode
            if retcode < 0:
                handle_error('app_terminated', result, mesh_name, stats_name)
            elif retcode > 0:
                handle_error('app_error', result, mesh_name, stats_name)
            else:
                result['ok'] += [mesh_name]
        except OSError:
            handle_error('os_error', result, mesh_name, stats_name)
        except subprocess.TimeoutExpired:
            handle_error('app_timeout', result, mesh_name, stats_name)
    with open(f'{dirname}/result.json', 'wt') as f:
        json.dump(result, f, indent=2)


@run.command()
@click.option('--dirname', '-d')
def draw(dirname):
    def handle_error(err, result, scene_name):
        msg = 'error: ' + err
        print(msg + ' ' * max(0, 78-len(msg)))
        if err not in result:
            result[err] = []
        result[err] += [scene_name]

    result = {'ok': []}
    scene_names = glob.glob(f'{dirname}/{scene_dir}/*.json')
    scene_num = len(scene_names)
    for scene_id, scene_name in enumerate(scene_names):
        image_name = scene_name.replace(
            'scenes/', 'images/').replace('.json', '.jpg')
        msg = f'[{scene_id}/{scene_num}] {scene_name}'
        print(msg + ' ' * max(0, 78-len(msg)))
        cmd = f'../yocto-gl/bin/yscenetrace {scene_name} -o {image_name} -t eyelight -r 640 -s 16'
        try:
            retcode = subprocess.run(cmd, timeout=5, shell=True).returncode
            if retcode < 0:
                handle_error('app_terminated', result, scene_name)
            elif retcode > 0:
                handle_error('app_error', result, scene_name)
            else:
                result['ok'] += [scene_name]
        except OSError:
            handle_error('os_error', result, scene_name)
        except subprocess.TimeoutExpired:
            handle_error('app_timeout', result, scene_name)
    with open(f'{dirname}/result.json', 'wt') as f:
        json.dump(result, f, indent=2)


@run.command()
@click.argument('dirname')
def render(dirname):
    pass


@run.command()
@click.argument('dirname')
def render():
    pass


if __name__ == '__main__':
    run()
