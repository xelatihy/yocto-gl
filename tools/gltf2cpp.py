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
// Parse int function.
inline void serialize_from_json(int& val, const json& js) {
    if (!js.is_number_integer()) throw runtime_error("integer expected");
    val = js;
}

// Parse float function.
inline void serialize_from_json(float& val, const json& js) {
    if (!js.is_number()) throw runtime_error("number expected");
    val = js;
}

// Parse bool function.
inline void serialize_from_json(bool& val, const json& js) {
    if (!js.is_boolean()) throw runtime_error("bool expected");
    val = js;
}

// Parse std::string function.
inline void serialize_from_json(string& val, const json& js) {
    if (!js.is_string()) throw runtime_error("string expected");
    val = js;
}

// Parse json function.
inline void serialize_from_json(json& val, const json& js) { val = js; }

// Parse support function.
template <typename T>
inline void serialize_from_json(T*& val, const json& js) {
    if (js.is_null()) { val = nullptr; return; }
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new T();
    serialize_from_json(*val, js);
}

// Parse support function.
template <typename T>
inline void serialize_from_json(vector<T>& vals, const json& js) {
    if (!js.is_array()) throw runtime_error("array expected");
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        serialize_from_json(v, js[i]);
        vals[i] = v;
    }
}

// Parse support function.
template <typename T>
inline void serialize_from_json(map<string, T>& vals, const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        serialize_from_json(vals[kv.key()], kv.value());
    }
}

// Parse support function.
template<typename T, size_t N>
inline void serialize_from_json(array<T, N>& vals, const json& js) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (N != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < N; i++) serialize_from_json(vals[i], js.at(i));
}

// Parse support function.
template<typename T, typename T1>
inline void serialize_from_json(T& val, const json& js, const vector<pair<T1, T>>& table) {
    auto v = T1();
    serialize_from_json(v, js);
    auto found = false;
    for(auto& kv : table) { if(kv.first == v) { val = kv.second; found = true; break; } }
    if(!found) throw runtime_error("bad enum value");
}

// Parse support function.
inline void serialize_from_json(vec2f& vals, const json& js) {
    serialize_from_json((array<float, 2>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(vec3f& vals, const json& js) {
    serialize_from_json((array<float, 3>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(vec4f& vals, const json& js) {
    serialize_from_json((array<float, 4>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(quat4f& vals, const json& js) {
    serialize_from_json((array<float, 4>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(mat4f& vals, const json& js) {
    serialize_from_json((array<float, 16>&)vals, js);
}

// Parse id function.
template <typename T>
inline void serialize_from_json(glTFid<T>& val, const json& js) {
    if (!js.is_number_integer()) throw runtime_error("int expected");
    val = glTFid<T>((int)js);
}

// Parses a glTFProperty object
inline void serialize_from_json(glTFProperty& val, const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
#if YGL_GLTFJSON
    if(js.count("extensions")) serialize_from_json(val.extensions, js.at("extensions"));
    if(js.count("extras")) serialize_from_json(val.extras, js.at("extras"));
#endif
}
'''

parse_fmt = '''
{{#types}}
{{#enums}}
// Parse a {{name}} enum
inline void serialize_from_json({{name}}& val, const json& js) {
    static vector<pair<{{item}}, {{name}}>> table = { {{#values}} { {{enum}}, {{name}}::{{label}} },{{/values}} };
    serialize_from_json(val, js, table);
}

{{/enums}}

// Parses a {{name}} object
inline void serialize_from_json({{name}}& val, const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
    {{#base}}serialize_from_json(({{base}}&)val, js);{{/base}}
    {{#properties}}{{^extension}}{{#required}}if (!js.count("{{name}}")) throw runtime_error("missing value");{{/required}}{{^required}}if (js.count("{{name}}")) {{/required}}serialize_from_json(val.{{name}}, js.at("{{name}}"));{{/extension}}{{/properties}}
    {{#has_extensions}}
    if (js.count("extensions")) {
        auto& js_ext = js["extensions"];
        {{#properties}}{{#extension}}if (js_ext.count("{{extension}}")) serialize_from_json(val.{{name}}, js_ext.at("{{extension}}"));{{/extension}}{{/properties}}
    }
    {{/has_extensions}}
}
{{/types}}
'''

dump_func = '''
// Converts int to json.
inline void serialize_to_json(int val, json& js) { js = val; }

// Converts float to json.
inline void serialize_to_json(float val, json& js) { js = val; }

// Converts bool to json.
inline void serialize_to_json(bool val, json& js) { js = val; }

// Converts string to json.
inline void serialize_to_json(const string& val, json& js) { js = val; }

// Converts json to json.
inline void serialize_to_json(const json& val, json& js) { js = val; }

// Dump support function.
template<typename T>
inline void serialize_to_json(const T* val, json& js) {
    if (!val) { js = nullptr; return; }
    if (!js.is_object()) js = json::object();
    serialize_to_json(*val, js);
}

// Dump support function.
template<typename T, size_t N>
inline void serialize_to_json(const array<T, N>& vals, json& js) {
    js = json::array();
    for (auto i = 0; i < N; i++) serialize_to_json(vals[i], js[i]);
}

// Dump support function.
template <typename T>
inline void serialize_to_json(const vector<T>& vals, json& js) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) serialize_to_json(vals[i], js[i]);
}

// Dump support function.
template <typename T>
inline void serialize_to_json(const map<string, T>& vals, json& js) {
    js = json::object();
    for (auto& kv : vals) serialize_to_json(kv.second, js[kv.first]);
}

// Dump support function.
template<typename T, typename T1>
inline void serialize_to_json(const T& val, json& js, const vector<pair<T1, T>>& table) {
    auto found = false;
    auto v = T1();
    for(auto& kv : table) { if(kv.second == val) { v = kv.first; found = true; break; } }
    if(!found) throw runtime_error("invalid value");
    serialize_to_json(v, js);
}

// Dump support function.
inline void serialize_to_json(const vec2f& vals, json& js) {
    serialize_to_json((const array<float, 2>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const vec3f& vals, json& js) {
    serialize_to_json((const array<float, 3>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const vec4f& vals, json& js) {
    serialize_to_json((const array<float, 4>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const quat4f& vals, json& js) {
    serialize_to_json((const array<float, 4>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const mat4f& vals, json& js) {
    serialize_to_json((const array<float, 16>&)vals, js);
}

// Converts glTFid to json.
template <typename T>
inline void serialize_to_json(const glTFid<T>& val, json& js) {
    js = (int)val;
}

// Check for default value
template <typename T>
inline bool operator==(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a == (int)b;
}

// Check for default value
template <typename T>
inline bool operator!=(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a != (int)b;
}

// Converts a glTFProperty object to JSON
inline void serialize_to_json(const glTFProperty& val, json& js) {
    if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
    if (!val.extensions.empty()) serialize_to_json(val.extensions, js["extensions"]);
    if (!val.extras.is_null()) dump_attr(val.extras, "extras", js);
#endif
}

'''

dump_fmt = '''
{{#types}}
{{#enums}}
// Converts a {{name}} enum to JSON
inline void serialize_to_json(const {{name}}& val, json& js) {
    static vector<pair<{{item}}, {{name}}>> table = { {{#values}} { {{enum}}, {{name}}::{{label}} },{{/values}} };
    serialize_to_json(val, js, table);
}

{{/enums}}

// Converts a {{name}} object to JSON
inline void serialize_to_json(const {{name}}& val, json& js) {
    static auto def = {{name}}();
    if (!js.is_object()) js = json::object();
    {{#base}}serialize_to_json((const {{base}}&)val, js);{{/base}}
    {{#properties}}{{^extension}}{{^required}}if (val.{{name}} != def.{{name}}) {{/required}}serialize_to_json(val.{{name}}, js["{{name}}"]);{{/extension}}{{/properties}}
    {{#properties}}{{#extension}}
    if ({{def_check}}) {
        auto& js_ext = js["extensions"];
        serialize_to_json(val.{{name}}, js_ext["{{extension}}"]);
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
            vjs['def_check'] = '!' + 'val.' + vjs['name'] + '.empty()'
        elif 'glTFid<' in vjs['type']:
            vjs['def_check'] = 'val.' + vjs['name'] + '.is_valid()'
        elif vjs['type'] in ['json']:
            vjs['def_check'] = '!' + 'val.' + vjs['name'] + '.is_null()'
        elif vjs['type'] in ['vec3f','vec2f','vec4f','quat4f','mat4f']:
            vjs['def_check'] = 'val.' + vjs['name'] + ' != ' + vjs['type'] + vjs['default']
        else:
            vjs['def_check'] = 'val.' + vjs['name'] + ' != ' + vjs['default']

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
