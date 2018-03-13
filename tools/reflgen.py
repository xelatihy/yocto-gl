#! /usr/bin/env python3 -b

# Generates type reflection data from C++ code and comments
# Works only in very limited cases

import json, os, pystache

def extract(cpp, key):
    ncpp = ''
    inside = False
    for line in cpp.splitlines(True):
        if line.startswith('// #codegen begin ' + key):
            inside = True
        elif line.startswith('// #codegen end ' + key):
            inside = False
        else:
            if inside: ncpp += line
    return ncpp

def substitute(cpp, val, key):
    ncpp = ''
    inside = False
    for line in cpp.splitlines(True):
        if line.startswith('// #codegen begin ' + key):
            ncpp += line
            ncpp += val
            inside = True
        elif line.startswith('// #codegen end ' + key):
            ncpp += line
            inside = False
        else:
            if not inside: ncpp += line
    return ncpp

def parse(cpp):
    types = []
    cur = None
    curm = None
    for line in cpp.splitlines():
        if line.isspace() or line.startswith('    // ') or line.startswith('        ') or line.startswith('    }'):
            continue
        elif line.startswith('    ~'):
            if not curm: continue
            if not curm['mem_name']:
                cur['members'].pop()
                curm = None
            continue
        elif line.startswith('/// '):
            if not cur:
                cur = { 'name': '', 'doc': '', 'enum': False, 'members': [] }
                types += [ cur ]
            cur['doc'] += line.replace('/// ', '').strip()+' '
        elif line.startswith('struct'):
            cur['name'] = line.replace('struct','').replace('{','').strip()
            cur['enum'] = False
        elif line.startswith('enum'):
            cur['name'] = line.replace('struct','').replace('enum','').replace('{','').strip()
            cur['enum'] = True
        elif line.startswith('};'):
            cur = None
        elif line.startswith('    /// '):
            if not curm:
                curm = {'mem_name': '', 'mem_doc': '', 'mem_def': '', 'mem_type': ''}
                cur['members'] += [ curm ]
            curm['mem_doc'] += line.replace('/// ', '').strip()+' '
        elif line.startswith('    ') and line[4] and curm and '~' not in line:
            if cur['enum']:
                curm['mem_name'] = line.replace(',','').strip().partition('=')[0].strip()
            else:
                curm['mem_name'] = line.replace(';','').strip().partition('=')[0].split()[-1].strip()
            curm = None
    enums = [ t for t in types if t['enum'] ]
    structs = [ t for t in types if not t['enum'] ]
    def parse_docval(doc,tag,default):
        if tag in doc:
            val = doc.partition(tag+'(')[2].partition(')')[0]
            doc = doc.partition(tag+'(')[0] + doc.partition(tag+'(')[2].partition(')')[2]
            return doc, val
        else:
            return doc, default
    for struct in structs:
        for mem in struct['members']:
            doc = mem['mem_doc']
            doc, mem['mem_semantic'] = parse_docval(doc,'@refl_semantic','value')
            doc, mem['mem_uilimits'] = parse_docval(doc,'@refl_uilimits','0,0')
            doc, mem['mem_shortname'] = parse_docval(doc,'@refl_shortname','')
            while '  ' in doc:
                doc = doc.replace('  ', ' ')
            doc = doc.strip()
            mem['mem_doc'] = doc
    return enums, structs

enum_names_fmt = '''
{{#enums}}
/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, {{name}}>>&
enum_names<{{name}}>() {
    static auto names = std::vector<std::pair<std::string, {{name}}>>{
        {{#members}}
        {"{{mem_name}}", {{name}}::{{mem_name}}},
        {{/members}}
    };
    return names;
}

{{/enums}}

{{#structs}}
/// Visit struct elements.
template <typename Visitor>
inline void visit({{name}}& val, Visitor&& visitor) {
    {{#members}}
    visitor(val.{{mem_name}}, visit_var{"{{mem_name}}",
        visit_var_type::{{mem_semantic}},"{{mem_doc}}",{{mem_uilimits}},"{{mem_shortname}}"});
    {{/members}}
}

{{/structs}}

'''

mustache = pystache.Renderer(escape=lambda x: x)

def gen_refl(cpp):
    enums, structs = parse(cpp)
    enum_names = mustache.render(enum_names_fmt, {'enums': enums}, {'structs': structs})
    return enum_names

filename = 'yocto/yocto_gl.h'

with open(filename) as f: cpp = f.read()

tags = []
for line in cpp.splitlines():
    if line.startswith('// #codegen begin refl-'):
        tags += [ line.replace('// #codegen begin refl-','').strip() ]
tags = list(set(tags))

for tag in tags:
    decl = extract(cpp, 'refl-'+tag)
    refl = gen_refl(decl)
    cpp = substitute(cpp, refl, 'reflgen-'+tag)
with open(filename, 'wt') as f: f.write(cpp)

os.system('/usr/local/bin/clang-format -i -style=file yocto/yocto_gl.h')
