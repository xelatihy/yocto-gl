#! /usr/bin/env python3 -b

import sys

infilename = sys.argv[1]
outfilename = sys.argv[2]

def rotate90(line):
    if line.startswith('v '):
        tokens = line.split()
        tokens[1], tokens[3] = tokens[3], tokens[1]
        return ' '.join(tokens) + '\n'
    if line.startswith('vn '):
        tokens = line.split()
        tokens[1], tokens[3] = tokens[3], tokens[1]
        return ' '.join(tokens) + '\n'
    return line

def rotate180(line):
    if line.startswith('v '):
        tokens = line.split()
        tokens[1] = ('-' + tokens[1]).replace('--','')
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    if line.startswith('vn '):
        tokens = line.split()
        tokens[1] = ('-' + tokens[1]).replace('--','')
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    return line

def rotate90yz(line):
    if line.startswith('v '):
        tokens = line.split()
        tokens[2], tokens[3] = tokens[3], tokens[2]
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    if line.startswith('vn '):
        tokens = line.split()
        tokens[2], tokens[3] = tokens[3], tokens[2]
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    return line

def rotate180yz(line):
    if line.startswith('v '):
        tokens = line.split()
        tokens[2] = ('-' + tokens[2]).replace('--','')
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    if line.startswith('vn '):
        tokens = line.split()
        tokens[2] = ('-' + tokens[2]).replace('--','')
        tokens[3] = ('-' + tokens[3]).replace('--','')
        return ' '.join(tokens) + '\n'
    return line

obj = ''
with open(infilename) as f:
    obj = f.read()

nobj = ''
for line in obj.splitlines(True):
    # nobj += rotate90(line)
    # nobj += rotate180(line)
    # nobj += rotate90yz(line)
    # nobj += rotate180(rotate180yz(line))
    nobj += rotate90yz(rotate180yz(line))

with open(outfilename, 'w') as f:
    f.write(nobj)
