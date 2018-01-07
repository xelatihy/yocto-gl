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
struct json_serialize_stack {
    vector<string> path = {"glTF"};
    string pathname() {
        auto p = std::string();
        for (auto n : path) p += '/' + n;
        return p;
    }
};

// Parse support function.
template <typename T>
inline void json_serialize_value(vector<T>& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        vals.resize(js.size());
        for (auto i = 0; i < js.size(); i++) {
            // this is contrived to support for vector<bool>
            auto v = T();
            json_serialize_value(v, js[i], reading, err);
            vals[i] = v;
        }
    } else {
        js = json::array();
        for (auto i = 0; i < vals.size(); i++) { json_serialize_value(vals[i], js[i], reading, err); }
    }
}

// Parse support function.
template <typename T>
inline void json_serialize_value(map<string, T>& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_object()) throw runtime_error("object expected");
        for (auto kv = js.begin(); kv != js.end(); ++kv) {
            json_serialize_value(vals[kv.key()], kv.value(), reading, err);
        }
    } else {
        js = json::object();
        for (auto&& kv : vals) { json_serialize_value(kv.second, js[kv.first], reading, err); }
    }
}

// Parses a pointer
template<typename T>
inline void json_serialize_value(T*& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_object()) throw runtime_error("object expected");
        if (!val) val = new T();
        json_serialize_value(*val, js, reading, err);
    } else {
        if (!js.is_object()) js = json::object();
        json_serialize_value(*val, js, reading, err);
    }
}

// Parse int function.
inline void json_serialize_value(int& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_number_integer()) throw runtime_error("integer expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse float function.
inline void json_serialize_value(float& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_number()) throw runtime_error("number expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse bool function.
inline void json_serialize_value(bool& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_boolean()) throw runtime_error("bool expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse std::string function.
inline void json_serialize_value(string& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_string()) throw runtime_error("string expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse json function.
inline void json_serialize_value(json& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        val = js;
    } else {
        js = val;
    }
}

// Parse support function.
inline void json_serialize_value(vec2f& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        if (2 != js.size()) throw runtime_error("wrong array size");
        for (auto i = 0; i < 2; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    } else {
        js = json::array();
        for (auto i = 0; i < 2; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    }
}

// Parse support function.
inline void json_serialize_value(vec3f& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        if (3 != js.size()) throw runtime_error("wrong array size");
        for (auto i = 0; i < 3; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    } else {
        js = json::array();
        for (auto i = 0; i < 3; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    }
}

// Parse support function.
inline void json_serialize_value(vec4f& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        if (4 != js.size()) throw runtime_error("wrong array size");
        for (auto i = 0; i < 4; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    } else {
        js = json::array();
        for (auto i = 0; i < 4; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    }
}

// Parse support function.
inline void json_serialize_value(quat4f& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        if (4 != js.size()) throw runtime_error("wrong array size");
        for (auto i = 0; i < 4; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    } else {
        js = json::array();
        for (auto i = 0; i < 4; i++) { json_serialize_value(vals[i], js[i], reading, err); }
    }
}

// Parse support function.
inline void json_serialize_value(mat4f& vals, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_array()) throw runtime_error("array expected");
        if (16 != js.size()) throw runtime_error("wrong array size");
        for (auto j = 0; j < 4; j++) {
            for (auto i = 0; i < 4; i++) { json_serialize_value(vals[j][i], js[j * 4 + i], reading, err); }
        }
    } else {
        js = json::array();
        for (auto j = 0; j < 4; j++) {
            for (auto i = 0; i < 4; i++) { json_serialize_value(vals[j][i], js[j * 4 + i], reading, err); }
        }
    }
}

// Parse support function.
template <typename T>
inline void json_serialize_attr(T& val, const char* name, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        auto iter = js.find(name);
        if (iter == js.end()) return;
        err.path.push_back(name);
        json_serialize_value(val, *iter, reading, err);
        err.path.pop_back();
    } else {
        err.path.push_back(name);
        json_serialize_value(val, js[name], reading, err);
        err.path.pop_back();
    }
}

// Parse support function.
template <typename T>
inline void json_serialize_enum(T& val, json& js, const vector<pair<string, T>>& vals, bool reading, json_serialize_stack& err) {
    if(reading) {
        auto v = string();
        json_serialize_value(v, js, reading, err);
        auto found = false;
        for(auto& kv : vals) { if(kv.first == v) { val = kv.second; found = true; break; } }
        if (!found) throw runtime_error("bad enum value");
    } else {
        auto v = string();
        for(auto& kv : vals) { if(kv.second == val) { v = kv.first; break; } }
        json_serialize_value(v, js, reading, err);
    }
}

// Parse support function.
template <typename T>
inline void json_serialize_enum(T& val, json& js, const vector<pair<int, T>>& vals, bool reading, json_serialize_stack& err) {
    if(reading) {
        auto v = 0;
        json_serialize_value(v, js, reading, err);
        auto found = false;
        for(auto& kv : vals) { if(kv.first == v) { val = kv.second; found = true; break; } }
        if (!found) throw runtime_error("bad enum value");
    } else {
        auto v = 0;
        for(auto& kv : vals) { if(kv.second == val) { v = kv.first; break; } }
        json_serialize_value(v, js, reading, err);
    }
}

// Parse id function.
template <typename T>
inline void json_serialize_value(glTFid<T>& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_number_integer()) throw runtime_error("int expected");
        val = glTFid<T>((int)js);
    } else {
        js = (int)val;
    }
}

// Parses a glTFProperty object
inline void json_serialize_value(glTFProperty& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_object()) throw runtime_error("object expected");
#if YGL_GLTFJSON
        json_serialize_attr(val.extensions, "extensions", js, err);
        json_serialize_attr(val.extras, "extras", js, err);
#endif
    } else {
        if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
        if (!val.extensions.empty())
            dump_attr(val.extensions, "extensions", js, err);
        if (!val.extras.is_null()) dump_attr(val.extras, "extras", js, err);
#endif
    }
}
'''

parse_fmt = '''
{{#types}}
{{#enums}}
// Parse a {{name}} enum
inline void json_serialize_value({{name}}& val, json& js, bool reading, json_serialize_stack& err) {
    static vector<pair<{{item}}, {{name}}>> table = { {{#values}} { {{enum}}, {{name}}::{{label}} },{{/values}} };
    json_serialize_enum(val, js, table, reading, err);
}

{{/enums}}

// Parses a {{name}} object
inline void json_serialize_value({{name}}& val, json& js, bool reading, json_serialize_stack& err) {
    if(reading) {
        if (!js.is_object()) throw runtime_error("object expected");
        {{#base}}json_serialize_value(({{base}}&)val, js, reading, err);{{/base}}
        {{#properties}}{{^extension}}{{#required}}if (!js.count("{{name}}")) throw runtime_error("missing required variable");{{/required}}json_serialize_attr(val.{{name}}, "{{name}}", js, reading, err);{{/extension}}{{/properties}}
        {{#has_extensions}}
        if (js.count("extensions")) {
            auto& js_ext = js["extensions"];
            {{#properties}}{{#extension}}json_serialize_attr(val.{{name}}, "{{extension}}", js_ext, reading, err);{{/extension}}{{/properties}}
        }
        {{/has_extensions}}
    } else {
        if (!js.is_object()) js = json::object();
        {{#base}}json_serialize_value(({{base}}&)val, js, reading, err);{{/base}}
        {{#properties}}{{^extension}}{{^required}}if ({{def_check}}) {{/required}}json_serialize_attr(val.{{name}}, "{{name}}", js, reading, err);{{/extension}}{{/properties}}
        {{#properties}}{{#extension}}
        if ({{def_check}}) {
            auto& js_ext = js["extensions"];
            json_serialize_attr(val.{{name}}, "{{extension}}", js_ext, reading, err);
        }
        {{/extension}}{{/properties}}
    }
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

# substitute
substitute('yocto/yocto_gl.h', types, 'type')
substitute('yocto/yocto_gl.cpp', funcs, 'func')
os.system('/usr/local/bin/clang-format -i -style=file yocto/yocto_gl.h')
os.system('/usr/local/bin/clang-format -i -style=file yocto/yocto_gl.cpp')
