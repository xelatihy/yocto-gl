#! /usr/bin/env python3 -B

import re

filename = 'bitterli/pbrt/cornell-box/scene.pbrt'

with open(filename) as f:
    pbrt = f.read()

pbrt = re.sub(r'"(\w+)\s+(\w+)"',r'"\1|\2"', pbrt)
pbrt = pbrt.replace('[',' [ ').replace(']', ' ] ')

last = ''
inparam = False
xpbrt = ''
next_name = False
for line in pbrt.splitlines():
    tokens = line.split()
    for tok in tokens:
        if tok == '[':
            inparam = True
            xpbrt += ' "'
        elif tok == ']':
            inparam = False
            xpbrt += '" '
        elif tok[0] == '"':
            if next_name:
                xpbrt += '"name"=' + tok + ' '
                next_name = False
            elif inparam:
                xpbrt += tok.replace('"','')
            else:
                xpbrt += tok + '='
        elif tok[0] in "0123456789+-":
            xpbrt += tok + ' '
        else:
            if last:
                xpbrt += '/>\n'
            elif tok == 'WorldBegin':
                xpbrt += '<World>\n'
            elif tok == 'WorldEnd':
                xpbrt += '</World>'
            elif tok == 'AttributeBegin':
                xpbrt += '<Attribute>\n'
            elif tok == 'AttributeEnd':
                xpbrt += '</Attribute>'
            else:
                last = tok
                xpbrt += '<' + tok + ' '
                next_name = tok not in ['Transform']
    xpbrt += '\n'

with open(filename+'.xml','wt') as f:
    f.write(xpbrt)

pbrt = pbrt.replace("Transform ", 'Transform "tranform" "float|xform" ')

last = ''
inparam = False
jpbrt = '[ \n'
next_name = False
for line in pbrt.splitlines():
    for tok in line.split():
        if tok == '[':
            inparam = True
            jpbrt += ' [ '
        elif tok == ']':
            inparam = False
            jpbrt += ' ], '
        elif tok[0] == '"':
            if next_name:
                jpbrt += '"name": ' + tok + ', '
                next_name = False
            elif inparam:
                jpbrt += tok + ', '
            else:
                jpbrt += tok + ': '
        elif tok[0] in "0123456789+-":
            jpbrt += tok + ', '
        else:
            if last:
                last = ''
                jpbrt += ' }, \n'
            if tok == 'WorldBegin':
                jpbrt += '{ "kind": "World", "content": [\n'
                continue
            elif tok == 'WorldEnd':
                jpbrt += '] }, '
                continue
            elif tok == 'AttributeBegin':
                jpbrt += '{ "kind": "Attribute", "content": [\n'
                continue
            elif tok == 'AttributeEnd':
                jpbrt += '] }, '
                continue
            else:
                last = tok
                jpbrt += '{ "kind": "' + tok + '", '
                next_name = True
jpbrt += '\n]\n\n'
jpbrt = jpbrt.replace(',  ]', ' ]').replace(',  }', ' }')

with open(filename+'.json','wt') as f:
    f.write(jpbrt)
