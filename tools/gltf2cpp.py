#! /usr/bin/env python3 -b

# Comvert types similar to JSON schema to C++ code

import json, os, pystache

type_fmt = '''
// forward declaration
{{#types}}
struct {{name}};
{{/types}}

{{#types}}

{{#enums}}
/// {{description}}
enum class {{name}} {
/// Not set
NotSet = -1,
{{#values}}
// {{description}}
{{label}} = {{value}},
{{/values}}
};

{{/enums}}

/// {{description}}
struct {{name}} {{#base}}: {{base}}{{/base}} {
{{#properties}}
/// {{description}}{{#required}} [required]{{/required}}
{{type}} {{name}} = {{default}};
{{/properties}}
{{#extra_properties}}
/// {{description}}{{#required}} [required]{{/required}}
{{type}} {{name}} = { };
{{/extra_properties}}
{{#has_getters}}{{#properties}}{{#is_ptrarray}}
/// typed access for nodes
{{item}}* get(const glTFid<{{item}}>& id) const {
    if (!id) return nullptr;
    return {{name}}.at((int)id);
}{{/is_ptrarray}}{{/properties}}{{/has_getters}}
{{#destructor}}
/// Cleanup
~{{name}}() {
{{#properties}}{{#is_ptr}}if ({{name}}) delete {{name}};{{/is_ptr}}{{#is_ptrarray}}for(auto v : {{name}}) if(v) delete v;{{/is_ptrarray}}{{/properties}}
}
{{/destructor}}
};
{{/types}}
'''

parse_func = '''
// Parse error
struct parse_stack {
    vector<string> path = {"glTF"};
    string pathname() {
        auto p = std::string();
        for (auto n : path) p += '/' + n;
        return p;
    }
};

// Parse support function.
template <typename T>
inline void parse(vector<T>& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        parse(v, js[i], err);
        vals[i] = v;
    }
}

// Parse int function.
inline void parse(int& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer()) throw runtime_error("integer expected");
    val = js;
}

// Parse float function.
inline void parse(float& val, const json& js, parse_stack& err) {
    if (!js.is_number()) throw runtime_error("number expected");
    val = js;
}

// Parse bool function.
inline void parse(bool& val, const json& js, parse_stack& err) {
    if (!js.is_boolean()) throw runtime_error("bool expected");
    val = js;
}

// Parse std::string function.
inline void parse(string& val, const json& js, parse_stack& err) {
    if (!js.is_string()) throw runtime_error("string expected");
    val = js;
}

// Parse json function.
inline void parse(json& val, const json& js, parse_stack& err) { val = js; }

// Parse support function.
inline void parse(vec2f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (2 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 2; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(vec3f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (3 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 3; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(vec4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (4 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 4; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(quat4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (4 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 4; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(mat4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (16 != js.size()) throw runtime_error("wrong array size");
    for (auto j = 0; j < 4; j++) {
        for (auto i = 0; i < 4; i++) { parse(vals[j][i], js[j * 4 + i], err); }
    }
}

// Parse support function.
template <typename T>
inline void parse(map<string, T>& vals, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        parse(vals[kv.key()], kv.value(), err);
    }
}

// Parse support function.
template <typename T>
inline void parse_attr(
    T& val, const char* name, const json& js, parse_stack& err) {
    auto iter = js.find(name);
    if (iter == js.end()) return;
    err.path.push_back(name);
    parse(val, *iter, err);
    err.path.pop_back();
}

// Parse id function.
template <typename T>
inline void parse(glTFid<T>& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer()) throw runtime_error("int expected");
    val = glTFid<T>((int)js);
}

// Parses a glTFProperty object
inline void parse(glTFProperty*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFProperty();
#if YGL_GLTFJSON
    parse_attr(val->extensions, "extensions", js, err);
    parse_attr(val->extras, "extras", js, err);
#endif
}
'''

parse_fmt = '''
{{#types}}
{{#enums}}
// Parse a {{name}} enum
inline void parse({{name}}& val, const json& js,
    parse_stack& err) {
    static map<{{item}}, {{name}}> table = { {{#values}} { {{enum}}, {{name}}::{{label}} },{{/values}} };
    auto v = {{item}}();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

{{/enums}}

// Parses a {{name}} object
inline void parse(
    {{name}}*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new {{name}}();
    {{#base}}parse(({{base}}*&)val, js, err);{{/base}}
    {{#properties}}{{^extension}}{{#required}}if (!js.count("{{name}}")) throw runtime_error("missing required variable");{{/required}}parse_attr(val->{{name}}, "{{name}}", js, err);{{/extension}}{{/properties}}
    {{#has_extensions}}
    if (js.count("extensions")) {
        auto& js_ext = js["extensions"];
        {{#properties}}{{#extension}}parse_attr(val->{{name}}, "{{extension}}", js_ext, err);{{/extension}}{{/properties}}
    }
    {{/has_extensions}}
}
{{/types}}
'''

dump_func = '''
// Dump support function.
template <typename T>
inline void dump(const vector<T>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) { dump(vals[i], js[i], err); }
}

// Converts int to json.
inline void dump(const int& val, json& js, parse_stack& err) { js = val; }

// Converts float to json.
inline void dump(const float& val, json& js, parse_stack& err) { js = val; }

// Converts bool to json.
inline void dump(const bool& val, json& js, parse_stack& err) { js = val; }

// Converts string to json.
inline void dump(const string& val, json& js, parse_stack& err) { js = val; }

// Converts json to json.
inline void dump(const json& val, json& js, parse_stack& err) { js = val; }

// Dump support function.
inline void dump(const vec2f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 2; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const vec3f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 3; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const vec4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 4; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const quat4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 4; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const mat4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto j = 0; j < 4; j++) {
        for (auto i = 0; i < 4; i++) { dump(vals[j][i], js[j * 4 + i], err); }
    }
}

// Dump support function.
template <typename T>
inline void dump(const map<string, T>& vals, json& js, parse_stack& err) {
    js = json::object();
    for (auto&& kv : vals) { dump(kv.second, js[kv.first], err); }
}

// Dump support function.
template <typename T>
inline void dump_attr(
    const T& val, const char* name, json& js, parse_stack& err) {
    err.path.push_back(name);
    dump(val, js[name], err);
    err.path.pop_back();
}

// Converts glTFid to json.
template <typename T>
inline void dump(const glTFid<T>& val, json& js, parse_stack& err) {
    js = (int)val;
}

// Converts a glTFProperty object to JSON
inline void dump(const glTFProperty* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
    if (!val->extensions.empty())
        dump_attr(val->extensions, "extensions", js, err);
    if (!val->extras.is_null()) dump_attr(val->extras, "extras", js, err);
#endif
}

'''

dump_fmt = '''
{{#types}}
{{#enums}}
// Converts a {{name}} enum to JSON
inline void dump(const {{name}}& val, json& js, parse_stack& err) {
    static map<{{name}}, {{item}}> table = { {{#values}} { {{name}}::{{label}}, {{enum}}  },  {{/values}} };
    auto v = table.at(val);
    dump(v, js, err);
}

{{/enums}}

// Converts a {{name}} object to JSON
inline void dump(const {{name}}* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    {{#base}}dump((const {{base}}*)val, js, err);{{/base}}
    {{#properties}}{{^extension}}{{^required}}if ({{def_check}}) {{/required}}dump_attr(val->{{name}}, "{{name}}", js, err);{{/extension}}{{/properties}}
    {{#properties}}{{#extension}}
    if ({{def_check}}) {
        auto& js_ext = js["extensions"];
        dump_attr(val->{{name}}, "{{extension}}", js_ext, err);
    }
    {{/extension}}{{/properties}}
}
{{/types}}
'''

def substitute(filename, val, key):
    with open(filename) as f: cpp = f.read()
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
    with open(filename, 'wt') as f: f.write(ncpp)

filenames = [
    "glTFChildOfRootProperty.schema.json",
    "accessor.sparse.indices.schema.json",
    "accessor.sparse.values.schema.json",
    "accessor.sparse.schema.json",
    "accessor.schema.json",
    "animation.channel.target.schema.json",
    "animation.channel.schema.json",
    "animation.sampler.schema.json",
    "animation.schema.json",
    "asset.schema.json",
    "buffer.schema.json",
    "bufferView.schema.json",
    "camera.orthographic.schema.json",
    "camera.perspective.schema.json",
    "camera.schema.json",
    "image.schema.json",
    "textureInfo.schema.json",
    "texture.schema.json",
    "material.normalTextureInfo.schema.json",
    "material.occlusionTextureInfo.schema.json",
    "material.pbrMetallicRoughness.schema.json",
    "material.pbrSpecularGlossiness.schema.json",
    "material.schema.json",
    "mesh.primitive.schema.json",
    "mesh.schema.json",
    "node.schema.json",
    "sampler.schema.json",
    "scene.schema.json",
    "skin.schema.json",
    "glTF.schema.json"
]

mustache = pystache.Renderer(escape=lambda x: x)

def fix_schema(js):
    print(js['name'])
    js['has_getters'] = js['name'] in ['glTF', 'glTFAnimation']
    for k, v in js['properties'].items(): v['name'] = k
    js['properties'] = [ v for _, v in js['properties'].items() if v ]
    if 'extra_properties' not in js: js['extra_properties'] = {}
    for k, v in js['extra_properties'].items(): v['name'] = k
    js['extra_properties'] = [ v for _, v in js['extra_properties'].items() if v ]
    if 'enums' not in js: js['enums'] = []
    js['destructor'] = False
    js['has_extensions'] = False
    for ejs in js['enums']:
        for vv in ejs['values']:
            if 'value' not in vv: vv['value'] = vv['enum']
            if isinstance(vv['enum'], str): vv['enum'] = '"' + vv['enum'] + '"'
            if 'description' not in vv: vv['description'] = vv['label']
    for vjs in js['properties']:
        vjs['required'] = 'required' in js and vjs['name'] in js['required']
        if '*' in vjs['type']:
            js['destructor'] = True
            if 'vector<' in vjs['type']: vjs['is_ptrarray'] = True
            else: vjs['is_ptr'] = True
        if 'extension' in vjs: js['has_extensions'] = True
        defaults = { 'int': '0', 'float': '0', 'string': '""' }
        if 'default' not in vjs:
            if 'vector' not in vjs['type'] and '*' in vjs['type']:
                vjs['default'] = 'nullptr'
            elif 'is_enum' in vjs:
                vjs['default'] = vjs['type']+'::NotSet'
            elif vjs['type'] in defaults:
                vjs['default'] = defaults[vjs['type']]
            else:
                vjs['default'] = '{}'
        elif 'is_enum' in vjs:
            vjs['default'] = vjs['type'] + '::' + vjs['default']
        else:
            vjs['default'] = str(vjs['default']).replace('.0','').replace('[','{').replace(']','}').replace('False','false')
        if 'vector<' in vjs['type'] or 'map<' in vjs['type']:
            vjs['def_check'] = '!' + 'val->' + vjs['name'] + '.empty()'
        elif 'glTFid<' in vjs['type']:
            vjs['def_check'] = 'val->' + vjs['name'] + '.is_valid()'
        elif vjs['type'] in ['json']:
            vjs['def_check'] = '!' + 'val->' + vjs['name'] + '.is_null()'
        elif vjs['type'] in ['vec3f','vec2f','vec4f','quat4f','mat4f']:
            vjs['def_check'] = 'val->' + vjs['name'] + ' != ' + vjs['type'] + vjs['default']
        else:
            vjs['def_check'] = 'val->' + vjs['name'] + ' != ' + vjs['default']

    return js

# load all JSONs
schemas = []
for filename in filenames:
    with open('tools/gltf/types/' + filename) as f:
        schemas += [ fix_schema(json.load(f)) ]

# make types
types = mustache.render(type_fmt, {'types': schemas})

# make funcs
funcs = ''
funcs += parse_func + '\n\n';
funcs += mustache.render(parse_fmt, {'types': schemas})
funcs += '\n\n';
funcs += dump_func + '\n\n';
funcs += mustache.render(dump_fmt, {'types': schemas})

# substitute
substitute('yocto/yocto_gl.h', types, 'type')
substitute('yocto/yocto_gl.cpp', funcs, 'func')
os.system('/usr/local/bin/clang-format -i -style=file yocto/yocto_gl.h')
os.system('/usr/local/bin/clang-format -i -style=file yocto/yocto_gl.cpp')
