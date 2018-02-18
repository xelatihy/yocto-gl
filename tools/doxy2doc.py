#! /usr/bin/env python3 -b

## This is not complete. Only here as reference.

import xml.etree.ElementTree as et
import json, subprocess, pystache, os, markdown
from collections import OrderedDict as odict

dir_xml = 'build/doxygen/xml'
dir_json = 'build/doxygen/json'
dir_md = 'build/doxygen/md'

def xml2dict(xml,
             array_tags=['memberdef','sectiondef','param'],
             text_tags=['briefdescription'], code_tags=['name','type','compoundname','argsstring','definition'],
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
    if tag in code_tags:
        code = (text + tail).strip()
        code = code.replace('< ','<').replace(' >','>').replace(' *','* ').replace(' &','& ').replace(',>','>')
        return tag, code
    return tag, (text + tail).strip()

def fix_template(js):
    if 'templateparamlist' not in js: return
    if not js['templateparamlist']:
        js['templatedecl'] = 'template<>'
    else:
        js['templatedecl'] = 'template<' + ', '.join(p['type'] + (p['declname'] if 'declname' in p else '') for p in js['templateparamlist']['param']) + '>'

def fix_type(js):
    if 'type' not in js: return
    if isinstance(js['type'], str):
        js['typestring'] = js['type']
    elif 'definition' in js:
        js['typestring'] = js['definition'].rpartition('ygl::')[0].strip()
    else:
        print('error')
        js['typestring'] = '<TYPE>'

def xml2json(xmlpath,jsonpath):
    xml = et.parse(xmlpath).getroot()
    _, js = xml2dict(xml)
    with open(jsonpath, 'wt') as f:
        json.dump(js, f, indent=2)

os.system('doxygen tools/Doxyfile')

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
            fix_template(member)
            fgroup[member['kind']] += [ member ]
    for struct in group['innerclasses']:
        fstruct = struct
        fix_template(fstruct)
        print(fstruct['compoundname'])
        if 'sectiondef' in fstruct:
            fstruct['function'] = [ ]
            fstruct['variable'] = [ ]
            for section in struct['sectiondef']:
                for member in section['memberdef']:
                    mname = member['name']
                    print(f'    {mname}')
                    print('        ' + member['kind'])
                    fix_template(member)
                    fix_type(member)
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
{{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}}
struct **{{compoundname}}**;

{{briefdescription}}

{{/declstruct}}

{{#struct}}
{{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}}
struct **{{compoundname}}**;

{{briefdescription}}

{{#function}}
- {{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}} {{type}} **{{name}}** {{argsstring}} <br>
  {{briefdescription}}

{{/function}}
{{#variable}}
- {{typestring}} **{{name}}** {{initializer}} <br>
  {{briefdescription}}
{{/variable}}

{{/struct}}

{{#enum}}
enum struct **{{name}}**;

{{briefdescription}}

{{#enumvalue}}
- **{{name}}** {{initializer}} <br>
  {{briefdescription}}
{{/enumvalue}}

{{/enum}}

{{#typedef}}
{{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}}
using **{{name}}** = {{type}}

{{briefdescription}}

{{/typedef}}

{{#variable}}
{{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}}
{{type}} **{{name}}** {{initializer}};

{{briefdescription}}

{{/variable}}

{{#function}}
{{#templatedecl}}{{templatedecl}} <br>{{/templatedecl}}
{{type}} **{{name}}**{{argsstring}};

{{briefdescription}}

{{/function}}
'''

api = '## API reference\n\n'
mustache = pystache.Renderer(escape=lambda x: x.replace('<','<').replace('>','\>').replace('*','\*').replace('_','\_'))
subprocess.run(f'rm -rf {dir_md} && mkdir -p {dir_md}', shell=True)
for fgroup in fgroups:
    name = fgroup['compoundname']
    md = mustache.render(doc_fmt, **fgroup)
    api += md
    with open(f'{dir_md}/{name}.md', 'wt') as f: f.write(md)

readme = ''
with open('yocto/yocto_gl.h') as f:
    for line in f:
        if line.startswith('/// '):
            readme += line[4:]
        elif line.startswith('///'):
            readme += line[3:]
        else:
            break
with open('readme.md', 'wt') as f: f.write(readme)

def add_toc(md):
    count = 0
    nmd = ''
    toc = ''
    for line in md.splitlines(True):
        if line.startswith('# '):
            title = 'About'
            toc += f'- [{title}](#toc{count})\n'
            nmd += f'<a id="toc{count}"></a>\n\n' + line
            count += 1
        elif line.startswith('## ') or line.startswith('### '):
            title = line.replace('#','').strip()
            indent = '    ' if line.startswith('### ') else ''
            toc += f'{indent}- [{title}](#toc{count})\n'
            nmd += f'<a id="toc{count}"></a>\n\n' + line
            count += 1
        else:
            nmd += line
    return nmd, toc

with open('docs/index.md','wt') as f:
    nmd, toc = add_toc(readme)
    f.write(nmd + toc)
with open('docs/api.md','wt') as f:
    nmd, toc = add_toc(readme + api)
    f.write(nmd + toc)

html_template = '''
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <title>Yocto/GL</title>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">
      <link rel="stylesheet" href="style.css">
    </head>
    <body>
    <header>
        <!-- a href="index.html">about</a -->
        <!-- a href="api.html">api</a -->
        <!-- a href="https://github.com/xelatihy/yocto-gl">github</a -->
        <p>
            <a href="index.html">Yocto/GL</a>
            <a href="https://github.com/xelatihy/yocto-gl">
            <img alt="" src="images/github-logo.png"></a>
            <a href="https://twitter.com/intent/tweet?text=Check%20out&amp;url=https%3A%2F%2Fgoo.gl%2FYvQvBr&amp;hashtags=yocto-gl&amp;via=xelatihy"><img alt="" src="images/twitter-logo.png"></a>
        </p>
        $toc$
    </header>
    <article>
    $body$
    <article>
    <footer></footer>
    </body>
    </html>
'''

def make_html(md):
    md, toc = add_toc(md)
    html = markdown.markdown(md, ['markdown.extensions.extra',
                     'markdown.extensions.codehilite'],
                     output_format='html5')
    htmltoc = markdown.markdown(toc, ['markdown.extensions.extra',
                     'markdown.extensions.codehilite'],
                     output_format='html5')
    html = html.replace('<pre>', '<pre><code>')
    html = html.replace('</pre>', '</code></pre>')
    # for link in glob.glob('docs/*.md'):
    #     link = link.replace('docs/','')
    #     hlink = link.replace('.md', '.html')
    #     html = html.replace(link, hlink)
    #     htmltoc = htmltoc.replace(link, hlink)
    while '<p><img' in html:
        before, _, remainder = html.partition('<p><img')
        middle, _, after = remainder.partition('</p>')
        html = before + '<figure><img' + middle + '</figure>' + after
    html = html_template.replace('$body$', html).replace('$toc$', htmltoc)
    return html

with open('docs/api.html','wt') as f:
    f.write(make_html(readme + api))
