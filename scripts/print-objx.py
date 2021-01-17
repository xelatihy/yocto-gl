#! /usr/bin/env python3 -b

import sys

filename = sys.argv[1]

with open(filename) as f:
    for line in f:
        if line.startswith('c'):
            tokens = line.split()
            print()
            print ('newmtl ycamera_{}'.format(tokens[1]))
            print ('  Ni {}'.format(float(tokens[3]) / float(tokens[4])))
            print ('  Ns {}'.format(round(float(tokens[5]) * 1000)))
            print ('  Kd {} {} {}'.format(float(tokens[8]), float(tokens[9]), float(tokens[10])))
            print ('  Ks {} {} {}'.format(float(tokens[8]) + float(tokens[14]) * float(tokens[6]), float(tokens[9]) + float(tokens[15]) * float(tokens[6]), float(tokens[10]) + float(tokens[16]) * float(tokens[6])))
        if line.startswith('e'):
            tokens = line.split()
            print()
            print ('newmtl yenvironment_{}'.format(tokens[1]))
            print ('  Ke {} {} {}'.format(tokens[2], tokens[3], tokens[4]))
            if tokens[5] != '""': print ('  map_Ke {}'.format(tokens[5]))
            print ('  Kd {} {} {}'.format(float(tokens[6]), float(tokens[7]), float(tokens[8])))
            print ('  Ks {} {} {}'.format(float(tokens[6]) + float(tokens[12]), float(tokens[7]) + float(tokens[13]), float(tokens[8]) + float(tokens[14])))
