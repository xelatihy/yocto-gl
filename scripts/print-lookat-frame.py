#! /usr/bin/env python3 -B

import glm, sys

if len(sys.argv) != 2 and len(sys.argv) != 3:
  print('usage: print-lookat-frame.py lookat [inv]')
  exit(1)

lookat_values = [ float(arg) for arg in sys.argv[1].replace('[', '').replace(']', '').split(',') ]

lookat_from = lookat_values[0:3]
lookat_to = lookat_values[3:6]
lookat_up = lookat_values[6:10]
invert = len(sys.argv) == 3

print('from:', lookat_from)
print('to: ', lookat_to)
print('up: ', lookat_up)
print('invert: ', invert)
