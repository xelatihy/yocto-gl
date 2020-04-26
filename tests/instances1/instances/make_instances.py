#! /usr/bin/env python3 -B

numx = 11
numy = 11

header = '''ply
format ascii 1.0
comment Written by Yocto/GL
comment https://github.com/xelatihy/yocto-gl
element instance {num}
property float xx
property float xy
property float xz
property float yx
property float yy
property float yz
property float zx
property float zy
property float zz
property float ox
property float oy
property float oz
end_header
'''

with open('instances.ply', 'w') as f:
  f.write(header.format(num=numx*numy))
  for j in range(numy):
    for i in range(numx):
      x = 0.2 * (i - (numx // 2))
      z = 0.2 * (j - (numy // 2))
      f.write(f'1 0 0 0 1 0 0 0 1 {x} 0 {z}\n')
