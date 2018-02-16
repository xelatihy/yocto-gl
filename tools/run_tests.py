#! /usr/bin/env python3 -b

import click, os, subprocess

test_dir = './tests'
out_dir = './run_tests/out'
check_dir = './run_tests/check'

qtrace_scenes = [
    'cornell_box.obj',
    'simple_al.obj',
    'lines_al.gltf',
    'points_al.gltf',
    'instancesl_pl.gltf',
]

def run_cmds(cmds):
    num = len(cmds)
    for idx, cmd in enumerate(cmds, 1):
        print(f'[CMD {idx}/{num}] {cmd}')
        ok = not subprocess.run(cmd, shell=True).returncode
        print(f'[OK  {idx}/{num}]' if ok else f'[ERR {idx}/{num}]')

@click.group()
def run():
    pass

@run.command()
@click.option('--copycheck', '-c', is_flag=True, default=False)
def qtrace(copycheck=False):
    odir = out_dir + '/qtrace'
    cdir = check_dir + '/qtrace'
    cmds = []
    cmds += [f'./tools/build.py release']
    cmds += [f'rm -rf {odir} && mkdir -p {odir}']
    if copycheck:
        cmds += [f'rm -rf {cdir} && mkdir -p {cdir}']        
    cmds += [f'./bin/ytestgen -q -c']
    for sname in qtrace_scenes:
        sdir = test_dir + '/' + sname.partition('.')[0]
        cmds += [ f'./bin/ytrace -q -s 16 -r 360 -o {odir}/{sname}.png {sdir}/{sname}' ]
        if copycheck:
            cmds += [ f'cp {odir}/{sname}.png {cdir}/{sname}.png' ]
        else:
            cmds += [ f'diff {odir}/{sname}.png {cdir}/{sname}.png' ]
    run_cmds(cmds)

if __name__ == '__main__':
    run()
