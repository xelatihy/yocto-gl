#! /usr/bin/env python3 -b

import click, os

test_dir = './tests'
out_dir = './tests/out'
check_dir = './tools/tests/out'

qtrace_scenes = [
    'simple_al.obj'
]

def run_cmds(cmds):
    num = len(cmds)
    for idx, cmd in enumerate(cmds):
        print(f'[{idx}/{num}] {cmd}')
        os.system(cmd)

@click.group()
def run():
    pass

@run.command()
def qtrace():
    odir = out_dir + '/qtrace'
    cdir = check_dir + '/qtrace'
    cmds = []
    cmds += [f'./tools/build.py release']
    cmds += [f'mkdir -p {odir}']
    for sname in qtrace_scenes:
        sdir = test_dir + '/' + sname.partition('.')[0]
        cmds += [ f'./bin/ytrace -s 16 -r 360 -o {odir}/{sname}.png {sdir}/{sname}' ]
    run_cmds(cmds)

if __name__ == '__main__':
    run()
