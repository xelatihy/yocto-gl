#! /usr/bin/env python3 -b

## This is not complete. Only here as reference.

import xml.etree.ElementTree as et
import textwrap

dirname = 'build/doxygen/xml/'

def unescape(txt):
    return txt.replace('&lt;','<').replace('&gt;','>').replace('&amp;','&').replace('&quot;','"').replace('&apos;','\'')

def descr2md(xml):
    txt = ''
    for t in xml.itertext():
        txt += t
    txt = unescape(txt)
    txt = txt.replace('\n',' ').replace('  ', ' ').strip()
    txt = textwrap.fill(txt)
    return txt

def doc2md(xml):
    return unescape(xml.text)

def code2md(text, indent='    '):
    text = unescape(text)
    text = text.replace(' *&','*& ').replace(' *','* ').replace(' &','& ').replace('< ','<').replace(' >','>').replace('* ,','*,')
    return ''.join(indent + l + '\n' for l in textwrap.wrap(text,78-len(indent)))

def template2md(xml):
    tmp_xml = elem_xml.find('templateparamlist')
    if tmp_xml:
        return code2md('template<' + ', '.join(e.text for e in tmp_xml.iter('type')) + '>')
    return ''

index_xml = et.parse(dirname + 'index.xml').getroot()
groups = [ e.attrib['refid'] for e in index_xml if e.attrib['kind'] == 'group' ]

md = ''
for group in groups:
    group_xml = et.parse(dirname + group + '.xml').getroot().find('compounddef')
    md += '### ' + doc2md(group_xml.find('title')) + '\n\n'
    for section_xml in group_xml.iter('sectiondef'):
        for elem_xml in section_xml.iter('memberdef'):
            if elem_xml.attrib['kind'] != 'enum': continue
            md += '#### Enum ' + doc2md(elem_xml.find('name')) + '\n\n'
            md += code2md('enum struct '+doc2md(elem_xml.find('name'))+' {')+'\n'
            for item_xml in elem_xml.iter('enumvalue'):
                if item_xml.find('initializer'):
                    md += code2md(item_xml.find('name').text + ' ' + item_xml.find('initializer').text + ',',indent='        ')
                else:
                    md += code2md(item_xml.find('name').text + ',', '        ')
            md += code2md('}')+'\n'
            md += descr2md(elem_xml.find('briefdescription')) + '\n\n'
            md += 'Members:\n'
            for item_xml in elem_xml.iter('enumvalue'):
                md += '    - ' + item_xml.find('name').text + ': ' + descr2md(item_xml.find('briefdescription')) +'\n'
            md += '\n'
    for class_ref in group_xml.iter('innerclass'):
        struct_xml = et.parse(dirname + class_ref.attrib['refid'] + '.xml').getroot()
        for elem_xml in section_xml.iter('compounddef'):
            md += '#### Struct ' + doc2md(elem_xml.find('compoundname')) + '\n\n'
            md += template2md(elem_xml.find('templateparamlist'))
            md += code2md('enum struct '+doc2md(elem_xml.find('name'))+' {')+'\n'
            for item_xml in elem_xml.iter('enumvalue'):
                if item_xml.find('initializer'):
                    md += code2md(item_xml.find('name').text + ' ' + item_xml.find('initializer').text + ',',indent='        ')
                else:
                    md += code2md(item_xml.find('name').text + ',', '        ')
            md += code2md('}')+'\n'
        pass
    for section_xml in group_xml.iter('sectiondef'):
        for elem_xml in section_xml.iter('memberdef'):
            if elem_xml.attrib['kind'] != 'typedef': continue
            md += '#### Typedef ' + doc2md(elem_xml.find('name')) + '\n\n'
            md += template2md(elem_xml.find('templateparamlist'))
            md += code2md(elem_xml.find('definition').text + ';').replace('typedef ', '') + '\n'
            md += descr2md(elem_xml.find('briefdescription')) + '\n\n'
    for section_xml in group_xml.iter('sectiondef'):
        for elem_xml in section_xml.iter('memberdef'):
            if elem_xml.attrib['kind'] != 'variable': continue
            md += '#### Constant ' + doc2md(elem_xml.find('name')) + '\n\n'
            md += code2md(elem_xml.find('definition').text + ' = ' + elem_xml.find('initializer').text + ';') + '\n'
            md += descr2md(elem_xml.find('briefdescription')) + '\n\n'
    for section_xml in group_xml.iter('sectiondef'):
        for elem_xml in section_xml.iter('memberdef'):
            if elem_xml.attrib['kind'] != 'function': continue
            md += '#### Function ' + doc2md(elem_xml.find('name')) + '()\n\n'
            md += template2md(elem_xml.find('templateparamlist'))
            md += code2md(elem_xml.find('definition').text + elem_xml.find('argsstring').text + ';') + '\n'
            md += descr2md(elem_xml.find('briefdescription')) + '\n\n'

with open('docs/doc.md', 'wt') as f: f.write(md)
