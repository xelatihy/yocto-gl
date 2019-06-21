#! /usr/bin/env python3 -B

from PIL import Image
import glob

for filename in glob.glob('*.png'):
    if '-thumb' in filename: continue
    img = Image.open(filename)
    aspect = img.size[0] / img.size[1]
    img = img.resize((740,int(740 / aspect)))
    img.save(filename.replace('.png', '-thumb.png'))
