#! /usr/bin/env python3 -B

import click

@click.group()
def run():
    pass

@run.command()
@click.argument('infilename')
@click.argument('outfilename')
def flipyz(infilename, outfilename):
    with open(infilename, 'rt') as f:
        obj = f.read()
    nobj = ''
    for line in obj.splitlines():
        if line.startswith('v ') or line.startswith('vn '):
            toks = line.split()
            toks[2], toks[3] = toks[3], toks[2]
            nobj += ' '.join(toks) + '\n'
        else:
            nobj += line + '\n'
    with open(outfilename, 'wt') as f:
        f.write(nobj)

@run.command()
@click.argument('infilename')
@click.argument('outfilename')
@click.option('--scale','-s',required=True)
def scale(infilename, outfilename, scale):
    with open(infilename, 'rt') as f:
        obj = f.read()
    nobj = ''
    for line in obj.splitlines():
        if line.startswith('v '):
            toks = line.split()
            toks[1] = str(float(toks[1]) * scale)
            toks[2] = str(float(toks[2]) * scale)
            toks[3] = str(float(toks[3]) * scale)
            nobj += ' '.join(toks) + '\n'
        else:
            nobj += line + '\n'
    with open(outfilename, 'wt') as f:
        f.write(nobj)

if __name__ == '__main__':
    run()
