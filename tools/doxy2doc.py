#! /usr/bin/env python3 -b

## This is not complete. Only here as reference.

import xml.etree.ElementTree as et
import json, subprocess, pystache
from collections import OrderedDict as odict

dir_xml = 'build/doxygen/xml'
dir_json = 'build/doxygen/json'
dir_md = 'build/doxygen/md'

def xml2dict(xml,
             array_tags=['memberdef','sectiondef','param'],
             text_tags=['briefdescription'],
             del_tags=['collaborationgraph','location','includes','detaileddescription','listofallmembers','inbodydescription']):
    tag = xml.tag
    text = xml.text if xml.text else ''
    tail = xml.tail if xml.tail else ''
    attr = xml.attrib
    children = list(xml)
    if tag in text_tags:
        s = text
        for st in xml.itertext():
            s += st
        s += tail
        s = s.replace('\n',' ').strip()
        return tag, s
    if children or attr:
        js = odict()
        for k, v in attr.items(): js[k] = v
        cjs = odict()
        for child in children:
            k, v = xml2dict(child)
            if k in del_tags: continue
            if k not in cjs: cjs[k] = []
            cjs[k] += [v]
        for k, v in cjs.items():
            if k in array_tags:
                js[k] = v
            else:
                js[k] = v if len(v) > 1 else v[0]
        return tag, js
    else:
        return tag, (text + tail).strip()

def xml2json(xmlpath,jsonpath):
    xml = et.parse(xmlpath).getroot()
    _, js = xml2dict(xml)
    with open(jsonpath, 'wt') as f:
        json.dump(js, f, indent=2)

index_xml = et.parse(f'{dir_xml}/index.xml').getroot()
group_names = [ e.attrib['refid'] for e in index_xml if e.attrib['kind'] == 'group' ]

subprocess.run(f'rm -rf {dir_json} && mkdir -p {dir_json}', shell=True)
groups = [ ]
for name in group_names:
    xml = et.parse(f'{dir_xml}/{name}.xml').getroot()
    group = xml2dict(xml)[1]['compounddef']
    if 'sectiondef' not in group:
        group['sectiondef'] = []
    class_names = [ e.attrib['refid'] for e in xml.iter('innerclass') ]
    if not class_names: continue
    group['innerclasses'] = [ ]
    for cname in class_names:
        cxml = et.parse(f'{dir_xml}/{cname}.xml').getroot()
        group['innerclasses'] += [ xml2dict(cxml)[1]['compounddef'] ]
    with open(f'{dir_json}/{name}.json', 'wt') as f: json.dump(group, f, indent=2)
    groups += [ group ]

fgroups = [ ]
for group in groups:
    name = group['compoundname']
    print(f'{name}')
    fgroup = odict()
    fgroup['compoundname'] = group['compoundname']
    fgroup['title'] = group['title']
    fgroup['enum'] = [ ]
    fgroup['declstruct'] = [ ]
    fgroup['struct'] = [ ]
    fgroup['function'] = [ ]
    fgroup['variable'] = [ ]
    fgroup['typedef'] = [ ]
    for section in group['sectiondef']:
        for member in section['memberdef']:
            mname = member['name']
            print(f'    {mname}')
            print('        ' + member['kind'])
            fgroup[member['kind']] += [ member ]
    for struct in group['innerclasses']:
        fstruct = struct
        print(fstruct['compoundname'])
        if 'sectiondef' in fstruct:
            fstruct['function'] = [ ]
            fstruct['variable'] = [ ]
            for section in struct['sectiondef']:
                for member in section['memberdef']:
                    mname = member['name']
                    print(f'    {mname}')
                    print('        ' + member['kind'])
                    fstruct[member['kind']] += [ member ]
            fgroup['struct'] += [ fstruct ]
            del fstruct['sectiondef']
        else:
            fgroup['declstruct'] += [ fstruct ]
    with open(f'{dir_json}/filtered__{name}.json', 'wt') as f: json.dump(fgroup, f, indent=2)
    fgroups += [ fgroup ]

doc_fmt = '''
### {{title}}

{{#declstruct}}
#### Struct {{compoundname}}

    struct {{name}}

{{briefdescription}}

{{/declstruct}}

{{#struct}}
#### Struct {{compoundname}}

    {{#templateparamlist}}template<{{#param}}{{type}}{{#declname}} {{declname}}{{/declname}}, {{/param}}>{{/templateparamlist}}
    struct {{name}} {
        {{#function}}
        {{#templateparamlist}}template<{{#param}}{{type}}{{#declname}} {{declname}}{{/declname}}, {{/param}}>{{/templateparamlist}}
        {{type}} {{name}}{{argsstring}};
        {{/function}}
        {{#variable}}
        {{name}}{{#initializer}} {{initializer}}{{/initializer}};
        {{/variable}}
    };

{{briefdescription}}

Members:
    {{#function}}
    - `{{name}}`: {{briefdescription}}
    {{/function}}
    {{#variable}}
    - `{{name}}`: {{briefdescription}}
    {{/variable}}

{{/struct}}

{{#enum}}
#### Enum `{{name}}`

    enum struct {{name}} {
        {{#enumvalue}}
        {{name}}{{#initializer}} {{initializer}}{{/initializer}},
        {{/enumvalue}}
    }

{{briefdescription}}

Members:
{{#enumvalue}}
    - `{{name}}`: {{briefdescription}}
{{/enumvalue}}

{{/enum}}

{{#typedef}}
#### Typedef `{{name}}`

    using {{name}} = {{type}}

{{briefdescription}}

{{/typedef}}

{{#variable}}
#### Constant `{{name}}`

    {{type}} {{name}} {{initializer}}

{{briefdescription}}

{{/variable}}

{{#function}}
#### Function `{{name}}()`

    {{#templateparamlist}}template<{{#param}}{{type}}{{#declname}} {{declname}}{{/declname}}, {{/param}}>{{/templateparamlist}}
    {{type}} {{name}}{{argsstring}};

{{briefdescription}}

{{/function}}
'''

def format_code(md):
    nmd = ''
    for line in md.splitlines():
        if line.startswith('    '):
            line = line.replace('< ','<').replace(' >','>').replace(' *','* ').replace(' &','& ').replace(',>','>')
        nmd += line + '\n'
    return nmd

mustache = pystache.Renderer(escape=lambda x: x)
subprocess.run(f'rm -rf {dir_md} && mkdir -p {dir_md}', shell=True)
for fgroup in fgroups:
    name = fgroup['compoundname']
    md = mustache.render(doc_fmt, **fgroup)
    md = format_code(md)
    with open(f'{dir_md}/{name}.md', 'wt') as f: f.write(md)
