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

def translatey(line,value):
    if line.startswith('v '):
        tokens = line.split()
        tokens[2] = str(float(tokens[2]) - value)
        return ' '.join(tokens) + '\n'
    return line

def scalexyz(line,value):
    if line.startswith('v '):
        tokens = line.split()
        tokens[1] = str(float(tokens[1]) * value)
        tokens[2] = str(float(tokens[2]) * value)
        tokens[3] = str(float(tokens[3]) * value)
        return ' '.join(tokens) + '\n'
    return line

def compute_bbox(line, bbox):
    if line.startswith('v '):
        tokens = line.split()
        vert = float(tokens[1]), float(tokens[2]), float(tokens[3])
        for index, coord in enumerate(vert):
            bbox[0][index] = min(bbox[0][index], vert[index])
            bbox[1][index] = max(bbox[1][index], vert[index])

obj = ''
with open(infilename) as f:
    obj = f.read()

inbbox = [ [1e10,1e10,1e10], [-1e10,-1e10,-1e10] ]
for line in obj.splitlines(True):
    compute_bbox(line, inbbox)
print(inbbox)

# exit()

nobj = ''
for line in obj.splitlines(True):
    # nobj += rotate90(line)
    # nobj += rotate180(line)
    # nobj += rotate90yz(line)
    # nobj += rotate180(rotate180yz(line))
    # nobj += rotate90yz(rotate180yz(line))
    # nobj += translatey(line, 3.0336)
    nobj += scalexyz(line, 0.1)
    # nobj += line

outbbox = [ [1e10,1e10,1e10], [-1e10,-1e10,-1e10] ]
for line in nobj.splitlines(True):
    compute_bbox(line, outbbox)
print(outbbox)

with open(outfilename, 'w') as f:
    f.write(nobj)
