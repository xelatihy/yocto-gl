///
/// YOCTO_GLTF_GEN: Code generation for Khronos GLTF loader and writer
/// form the schema.
///
/// This code is realsed only as a demonstration for others. It should not
/// be considered stable or even a good model for coding standards.
/// Unfortunately glTF changes too much and too often to have a stable codegen
/// path.
///

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../yocto/ext/json.hpp"
using nlohmann::json;

#include "../yocto/yocto_utils.h"

using namespace yu::operators;
using namespace yu::string;
using namespace yu::file;
using namespace yu::cmdline;

std::string filename2schemaname(const std::string& filename) {
    auto name = filename;
    name = replace(name, ".json", "");
    name = replace(name, ".schema", "");
    return name;
}

std::pair<std::string, std::string> split_schemaname(
    const std::string& schemaname) {
    auto pos = schemaname.rfind(".");
    if (pos == schemaname.npos) return {schemaname, ""};
    return {schemaname.substr(pos + 1), schemaname.substr(0, pos)};
}

std::string schemaname2typename(const std::string& schemaname) {
    if (startswith(schemaname, "glTF")) return schemaname;
    auto sname = schemaname;
    auto name = std::string();
    while (!sname.empty()) {
        auto split = split_schemaname(sname);
        name = (char)toupper(split.first[0]) + split.first.substr(1) + name;
        sname = split.second;
    }
    return "glTF" + name;
}

std::string varname2enumname(const std::string& varname) {
    return (char)toupper(varname[0]) + varname.substr(1);
}

std::string string2enumvalue(const std::string& str) {
    auto name = replace(lower(str), "/", "_");
    auto nname = std::string();
    while (!name.empty()) {
        auto pos = name.find("_");
        nname += (char)toupper(name[0]) + name.substr(1, pos - 1);
        name = (pos == name.npos) ? "" : name.substr(pos + 1);
    }
    name = nname;
    return name;
}

enum ValidCheckType { def, min, max, minLength, maxLength };

struct Type {
    virtual std::string sname() { return ""; }
    virtual std::string tname() { return ""; }
    virtual std::string fname() { return tname(); }
    virtual std::string rname() { return tname(); }
    virtual std::string to_type() { return ""; }
    virtual std::string to_default(const json& js) { return "{}"; }
    virtual std::string to_parse_func() { return ""; }
    virtual std::string to_dump_func() { return ""; }
    virtual std::string to_equal_func() { return ""; }
    virtual std::string to_validate_func() { return ""; }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "";
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return "";
    }

    virtual bool is_id() { return false; }
    virtual bool is_struct() { return false; }
    virtual bool is_array() { return false; }

    static std::string to_parse_support_func() {
        return R"___(
        // Parse error
        struct parse_stack {
            std::vector<std::string> path = {"glTF"};
            std::string pathname() {
                auto p = std::string();
                for(auto n : path) p += '/' + n;
                return p;
            }
        };

        // Parse support function.
        template <typename T>
        static bool parse(std::vector<T>& vals, const json& js,
                                  parse_stack& err) {
            if (!js.is_array()) return false;
            vals.resize(js.size());
            for (auto i = 0; i < js.size(); i++) {
                // this is contrived to support for vector<bool>
                auto v = T();
                if(!parse(v, js[i], err)) return false;
                vals[i] = v;
            }
            return true;
        }

        // Parse support function.
        template <typename T, int N>
        static bool parse(ym::vec<T, N>& vals, const json& js,
                                  parse_stack& err) {
            if (!js.is_array()) return false;
            if (N != js.size()) return false;
            for (auto i = 0; i < N; i++) {
                if(!parse(vals[i], js[i], err)) return false;
            }
            return true;
        }

        // Parse support function.
        template <typename T, int N>
        static bool parse(ym::quat<T, N>& vals, const json& js,
                                  parse_stack& err) {
            if (!js.is_array()) return false;
            if (N != js.size()) return false;
            for (auto i = 0; i < N; i++) {
                if(!parse(vals[i], js[i], err)) return false;
            }
            return true;
        }

        // Parse support function.
        template <typename T, int N, int M>
        static bool parse(ym::mat<T, N, M>& vals, const json& js,
                                  parse_stack& err) {
            if (!js.is_array()) return false;
            if (N*M != js.size()) return false;
            for (auto j = 0; j < M; j++) {
                for (auto i = 0; i < N; i++) {
                    if(!parse(vals[j][i], js[j*N+i], err)) return false;
                }
            }
            return true;
        }

        // Parse support function.
        template <typename T>
        static bool parse(std::map<std::string, T>& vals, const json& js,
                                  parse_stack& err) {
            if (!js.is_object()) return false;
            for (auto kv = js.begin(); kv != js.end(); ++kv) {
                if(!parse(vals[kv.key()], kv.value(), err)) return false;
            }
            return true;
        }

        // Parse support function.
        template <typename T>
        static bool parse_attr(T& val, const char* name,
                               const json& js, parse_stack& err) {
            auto iter = js.find(name);
            if(iter == js.end()) return true;
            err.path.push_back(name);
            if(!parse(val, *iter, err)) return false;
            err.path.pop_back();
            return true;
        }

        )___";
    }

    static std::string to_dump_support_func() {
        return R"___(
            // Dump support function.
            template <typename T>
            static void dump(const std::vector<T>& vals, json& js,
                                     parse_stack& err) {
                js = json::array();
                for (auto i = 0; i < vals.size(); i++) {
                    dump(vals[i], js[i], err);
                }
            }

            // Dump support function.
            template <typename T, int N>
            static void dump(const ym::vec<T, N>& vals, json& js,
                                     parse_stack& err) {
                js = json::array();
                for (auto i = 0; i < N; i++) {
                    dump(vals[i], js[i], err);
                }
            }

            // Dump support function.
            template <typename T, int N>
            static void dump(const ym::quat<T, N>& vals, json& js,
                                     parse_stack& err) {
                js = json::array();
                for (auto i = 0; i < N; i++) {
                    dump(vals[i], js[i], err);
                }
            }

            // Dump support function.
            template <typename T, int N, int M>
            static void dump(const ym::mat<T, N, M>& vals, json& js,
                                     parse_stack& err) {
                js = json::array();
                for (auto j = 0; j < M; j++) {
                    for (auto i = 0; i < N; i++) {
                        dump(vals[j][i], js[j*N+i], err);
                    }
                }
            }

            // Dump support function.
            template <typename T>
            static void dump(const std::map<std::string, T>& vals, json& js,
                                      parse_stack& err) {
                js = json::object();
                for (auto&& kv : vals) {
                    dump(kv.second, js[kv.first], err);
                }
            }

            // Dump support function.
            template <typename T>
            static void dump_attr(const T& val, const char* name,
                                              json& js, parse_stack& err) {
                err.path.push_back(name);
                dump(val, js[name], err);
                err.path.pop_back();
            }

        )___";
    }

    static std::string to_validate_support_func() {
        return R"___(
        // Validate support function.
        template <typename T>
        static void validate(const std::vector<T>& vals, parse_stack& err,
                             std::vector<std::pair<std::string, std::string>>& errs) {
            for (auto i = 0; i < vals.size(); i++) {
                validate(vals[i], err, errs);
            }
        }

        // Validate support function.
        template <typename T, int N>
        static void validate(const ym::vec<T, N>& vals, parse_stack& err,
                         std::vector<std::pair<std::string, std::string>>& errs) {
        }
        
        // Validate support function.
        template <typename T, int N>
        static void validate(const ym::quat<T, N>& vals, parse_stack& err,
                         std::vector<std::pair<std::string, std::string>>& errs) {
        }
        
        // Validate support function.
        template <typename T, int N, int M>
        static void validate(const ym::mat<T, N, M>& vals, parse_stack& err,
                         std::vector<std::pair<std::string, std::string>>& errs) {
        }
        
        // Validate support function.
        template <typename T>
        static void validate(const std::map<std::string, T>& vals, parse_stack& err,
                         std::vector<std::pair<std::string, std::string>>& errs) {
            for (auto&& kv : vals) {
                if(kv.first == "") errs.push_back({"missing dictionary key", err.pathname()});
                validate(kv.second, err, errs);
            }
        }

        // Validate support function.
        template <typename T>
        static void validate_attr(const T& val, const char* name, parse_stack& err,
                              std::vector<std::pair<std::string, std::string>>& errs) {
            err.path.push_back(name);
            validate(val, err, errs);
            err.path.pop_back();
        }

        )___";
    }
};

struct ExternType : Type {
    std::string schema_name = "";
    std::string type_name = "";

    ExternType(const std::string& sname, const std::string& tname)
        : schema_name(sname), type_name(tname) {}

    virtual std::string tname() { return type_name; }
    virtual std::string sname() { return schema_name; }
};

struct BuiltinType : Type {
    virtual std::string to_dump_func() {
        auto tmp = R"___(
            // Converts __TYPE__ to json.
            static void dump(const __TYPE__& val, json& js, parse_stack& err) {
                js = val;
            }
        )___";
        return replace(tmp, "__TYPE__", tname());
    }
    virtual std::string to_validate_func() {
        auto tmp = R"___(
        // Validates __TYPE__ (placeholder).
        static void validate(const __TYPE__& val, parse_stack& err,
                             std::vector<std::pair<std::string, std::string>>& errs) {
        }
        )___";
        return replace(tmp, "__TYPE__", tname());
    }
};

struct IntType : BuiltinType {
    virtual std::string tname() { return "int"; }
    virtual std::string sname() { return "integer"; }
    virtual std::string to_parse_func() {
        return R"___(
            // Parse int function.
            static bool parse(int& val, const json& js, parse_stack& err) {
                if (!js.is_number_integer()) return false;
                val = js;
                return true;
            }
        )___";
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) return "0";
        return std::to_string(js.get<int>());
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname + " != " + to_default(def);
    }
};

struct IdType : BuiltinType {
    virtual bool is_id() { return true; }
    virtual std::string tname() { return "glTFid<void>"; }
    virtual std::string sname() { return filename2schemaname("glTFid"); }
    virtual std::string to_parse_func() {
        return R"___(
        // Parse id function.
        template<typename T>
        static bool parse(glTFid<T>& val, const json& js, parse_stack& err) {
            if (!js.is_number_integer()) return false;
            val = glTFid<T>((int)js);
            return true;
        }
        )___";
    }
    virtual std::string to_dump_func() {
        return R"___(
        // Converts __TYPE__ to json.
        template<typename T>
        static void dump(const glTFid<T>& val, json& js, parse_stack& err) {
            js = (int)val;
        }
        )___";
    }
    virtual std::string to_validate_func() {
        return R"___(
        // Converts __TYPE__ to json.
        template<typename T>
        static void validate(const glTFid<T>& val, parse_stack& err,
                             std::vector<std::pair<std::string,std::string>>& errs) {
        }
        )___";
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) return "{}";
        return std::to_string(js.get<int>());
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "(bool)" + varname;
    }
};

struct TypedIdType : BuiltinType {
    std::string name = "";

    TypedIdType(const std::string& name) : name(name) {}

    virtual std::string tname() { return "glTFid<" + name + ">"; }
    virtual std::string sname() { return ""; }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "(bool)" + varname;
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return "!(bool)" + varname;
    }
};

struct StringType : BuiltinType {
    virtual std::string tname() { return "std::string"; }
    virtual std::string sname() { return "string"; }
    virtual std::string to_invalid_check() { return " == \"\""; }
    virtual std::string to_parse_func() {
        return R"___(
            // Parse std::string function.
            static bool parse(std::string& val, const json& js, parse_stack& err) {
                if (!js.is_string())  return false;
                val = js;
                return true;
            }
        )___";
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) return "\"\"";
        return "\"" + js.get<std::string>() + "\"";
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "!" + varname + ".empty()";
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return varname + ".empty()";
    }
};

struct FloatType : BuiltinType {
    virtual std::string tname() { return "float"; }
    virtual std::string sname() { return "number"; }
    virtual std::string to_parse_func() {
        return R"___(
            // Parse float function.
            static bool parse(float& val, const json& js, parse_stack& err) {
                if (!js.is_number())  return false;
                val = js;
                return true;
            }
        )___";
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) return "0";
        char buf[4096];
        sprintf(buf, "%g", js.get<float>());
        return buf;
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname + " != " + to_default(def);
    }
};

struct BoolType : BuiltinType {
    virtual std::string tname() { return "bool"; }
    virtual std::string sname() { return "boolean"; }
    virtual std::string to_parse_func() {
        return R"___(
            // Parse bool function.
            static bool parse(bool& val, const json& js, parse_stack& err) {
                if (!js.is_boolean())  return false;
                val = js;
                return true;
            }
        )___";
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) return "false";
        return (js.get<bool>()) ? "true" : "false";
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname;
    }
};

struct JsonType : BuiltinType {
    virtual std::string tname() { return "json"; }
    virtual std::string sname() { return "object"; }
    virtual std::string to_parse_func() {
        return R"___(
            // Parse json function.
            static bool parse(json& val, const json& js, parse_stack& err) {
                val = js;
                return true;
            }
        )___";
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "!" + varname + ".empty()";
    }
};

struct AliasType : Type {
    std::string schema_name = "";
    std::string type_name = "";
    Type* alias = nullptr;

    AliasType(const std::string& sname, const std::string& tname, Type* alias)
        : schema_name(sname), type_name(tname), alias(alias) {}

    virtual std::string tname() { return type_name; }
    virtual std::string sname() { return schema_name; }
    virtual std::string to_default(const json& js) {
        return alias->to_default(js);
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return alias->to_defined_value(varname, def);
    }
};

struct ArrayType : Type {
    Type* items = nullptr;

    ArrayType(Type* items) : items(items) {}

    virtual bool is_array() { return true; }

    virtual std::string tname() {
        return "std::vector<" + items->rname() + ">";
    }
    virtual std::string to_invalid_check() { return ".empty()"; }
    virtual std::string to_default(const json& js) {
        auto cpp = std::string();
        cpp += "{";
        if (!js.empty()) {
            auto count = 0;
            for (auto val_it = js.begin(); val_it != js.end(); ++val_it) {
                if (count++) cpp += ", ";
                cpp += items->to_default(*val_it);
            }
        }
        cpp += "}";
        return cpp;
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "!" + varname + ".empty()";
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return varname + ".empty()";
    }
    virtual std::string to_check(ValidCheckType type, const json& js) {
        if (type == ValidCheckType::def) return "";
        if (js.empty()) return "";
        switch (type) {
            case ValidCheckType::minLength:
                return ".size() >= " + std::to_string(js.get<int>());
            case ValidCheckType::maxLength:
                return ".size() <= " + std::to_string(js.get<int>());
            default: return "";
        };
    }
};

struct FixedArrayType : Type {
    std::string name;
    int nitems1 = 0;
    int nitems2 = 0;

    FixedArrayType(const std::string& name, int nitems1, int nitems2)
        : name(name), nitems1(nitems1), nitems2(nitems2) {}

    virtual std::string tname() { return name; }
    virtual std::string to_default(const json& js) {
        auto cpp = std::string();
        cpp += "{";
        auto defaults = std::vector<std::string>();
        if (js.empty()) {
            for (auto i = 0; i < nitems1; i++) {
                if (nitems2) {
                    for (auto j = 0; j < nitems2; j++) defaults.push_back("0");
                } else {
                    defaults.push_back("0");
                }
            }
        } else {
            auto items = FloatType();
            for (auto val_it = js.begin(); val_it != js.end(); ++val_it) {
                defaults.push_back(items.to_default(*val_it));
            }
        }
        auto count1 = 0;
        for (auto i = 0; i < nitems1; i++) {
            if (count1++) cpp += ", ";
            if (nitems2) {
                auto count2 = 0;
                cpp += "{";
                for (auto j = 0; j < nitems2; j++) {
                    if (count2++) cpp += ", ";
                    cpp += defaults[i * nitems2 + j];
                }
                cpp += "}";
            } else {
                cpp += defaults[i];
            }
        }
        cpp += "}";
        return cpp;
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname + " != " + tname() + to_default(def);
    }
};

struct Variable {
    std::string name = "";
    std::string help = "";
    std::string long_help = "";
    Type* type = nullptr;
    bool required = false;
    json defvalue = json();
    std::string ext_name = "";

    Variable(const std::string& name, Type* type) : name(name), type(type) {}

    virtual std::string to_type() {
        auto cpp = std::string();
        cpp += "    /// ";
        cpp += (help.empty()) ? "No documentation in schema." : help;
        if (required) cpp += " [required] ";
        cpp += "\n";
        //        if(long_help != "") cpp += "///\n/// " +
        //        replace_str(long_help, "\n", "\n///") + "\n";
        cpp += "    " + type->rname() + " " + name + " = " +
               type->to_default(defvalue) + ";\n";
        return cpp;
    }
};

struct StructType : Type {
    std::string schema_name = "";
    std::string type_name = "";
    std::string help = "";
    std::vector<Variable*> vars = {};
    std::vector<Variable*> extra_vars = {};
    StructType* base = nullptr;
    std::vector<Type*> subtypes;

    StructType(const std::string& sname, const std::string& tname)
        : schema_name(sname), type_name(tname) {}

    virtual bool is_struct() { return true; }
    virtual bool is_child_of_root() {
        if (!base) return false;
        return base->sname() == "glTFChildOfRootProperty";
    }

    virtual std::string tname() { return type_name; }
    virtual std::string sname() { return schema_name; }
    virtual std::string rname() { return tname() + "*"; }
    virtual std::string to_type() {
        auto cpp = std::string();
        cpp += "///\n/// " +
               ((help.empty()) ? "No documentation in schema." : help) +
               "\n///\n";
        cpp += "struct " + tname();
        if (base) cpp += " : " + base->tname();
        cpp += " {\n";
        for (auto var : vars) cpp += var->to_type();
        if (!extra_vars.empty()) {
            cpp += "\n// extra vars -------------------\n";
            for (auto var : extra_vars) cpp += var->to_type();
        }
        if (sname() == "glTF" || sname() == "animation") {
            cpp += "\n// access functions -------------\n";
            // cpp += "    // clang-format off\n";
            for (auto var : vars) {
                if (!var->type->is_array()) continue;
                auto atype = (ArrayType*)var->type;
                auto itype = atype->items;
                if (!itype->is_struct()) continue;
                cpp += "/// typed access for " + var->name + "\n";
                cpp += itype->tname() + "* get(const glTFid<" + itype->tname() +
                       ">& id) const {\n";
                cpp += "    if(!id) return nullptr;\n";
                cpp += "    return " + var->name + ".at( (int)id );\n";
                cpp += "}\n";
            }
            // cpp += "    // clang-format on\n";
        }
        cpp += to_destructor_func();
        cpp += "};\n\n";
        return cpp;
    }
    virtual std::string to_parse_func() {
        auto cpp = std::string();
        cpp += "// Parses a " + fname() + " object\n";
        cpp += "static bool parse(" + fname() +
               "*& val, const json& js, parse_stack& err) {\n";
        cpp += "    if(!js.is_object()) return false;\n";
        cpp += "    if(!val) val = new " + fname() + ";\n";
        if (base) {
            cpp += "    if(!parse((" + base->fname() +
                   "*&)val, js, err)) return false;\n";
        }
        for (auto var : vars) {
            if (var->ext_name != "") continue;
            if (var->required) {
                cpp += "    if(!js.count(\"" + var->name + "\")) return false;";
                cpp += "    if(!parse_attr(val->" + var->name + ", \"" +
                       var->name + "\", js, err)) return false;\n";
            } else {
                cpp += "    if(!parse_attr(val->" + var->name + ", \"" +
                       var->name + "\", js, err)) return false;\n";
            }
        }
        auto has_ext = false;
        for (auto var : vars)
            if (var->ext_name != "") has_ext = true;
        if (has_ext) {
            cpp += "if(js.count(\"extensions\")) {\n";
            cpp += "    auto& js_ext = js[\"extensions\"];";
            for (auto var : vars) {
                if (var->ext_name == "") continue;
                cpp += "    parse_attr(val->" + var->name + ", \"" +
                       var->ext_name + "\", js_ext, err);\n";
            }
            cpp += "}\n";
        }
        cpp += " return true;\n";
        cpp += "}\n\n";
        return cpp;
    }

    virtual std::string to_dump_func() {
        auto cpp = std::string();
        cpp += "// Converts a " + tname() + " object to JSON\n";
        cpp += "static void dump(const " + tname() +
               "* val, json& js, parse_stack& err) {\n";
        cpp += "    if (!js.is_object()) js = json::object();\n";
        if (base) {
            cpp += "    dump((const " + base->tname() + "*)val, js, err);\n";
        }
        for (auto var : vars) {
            if (var->ext_name != "") continue;
            if (var->required) {
                cpp += "    dump_attr(val->" + var->name + ", \"" + var->name +
                       "\", js, err);\n";
            } else {
                cpp += "    if(" +
                       var->type->to_defined_value(
                           "val->" + var->name, var->defvalue) +
                       ")";
                cpp += "    dump_attr(val->" + var->name + ", \"" + var->name +
                       "\", js, err);\n";
            }
        }
        auto has_ext = false;
        for (auto var : vars)
            if (var->ext_name != "") has_ext = true;
        if (has_ext) {
            cpp += "if(";
            auto count = 0;
            for (auto var : vars) {
                if (var->ext_name == "") continue;
                if (count++) cpp += " || ";
                cpp += var->type->to_defined_value(
                    "val->" + var->name, var->defvalue);
            }
            cpp += ") {";
            cpp += "    auto& js_ext = js[\"extensions\"];";
            for (auto var : vars) {
                if (var->ext_name == "") continue;
                cpp += "    if(" +
                       var->type->to_defined_value(
                           "val->" + var->name, var->defvalue) +
                       ")";
                cpp += "    dump_attr(val->" + var->name + ", \"" +
                       var->ext_name + "\", js_ext, err);\n";
            }
            cpp += "}";
        }
        cpp += "}\n\n";
        return cpp;
    }
    virtual std::string to_validate_func() {
        auto cpp = std::string();
        cpp += "// Validates a " + tname() + " object\n";
        cpp += "static void validate(const " + tname() +
               "* val, parse_stack& err, "
               "std::vector<std::pair<std::string,std::string>>& errs) {\n";
        cpp += "     if(!val) return;";
        if (base) {
            cpp +=
                "    validate((const " + base->tname() + "*)val, err, errs);\n";
        }
        for (auto var : vars) {
            if (var->required) {
                auto check = var->type->to_invalid_value(
                    "val->" + var->name, var->defvalue);
                if (!check.empty())
                    cpp += "    if(" + check +
                           ") errs.push_back({\"missing requried value " +
                           var->name + "\",err.pathname()});\n";
            }
            cpp += "    validate_attr(val->" + var->name + ", \"" + var->name +
                   "\", err, errs);\n";
        }
        cpp += "}\n\n";
        return cpp;
    }
    virtual std::string to_equal_func() {
        auto cpp = std::string();
        cpp +=
            "// Equality check for defaults (might go away in the "
            "future)\n";
        cpp += "static bool operator==(const " + tname() + "& a, const " +
               tname() + "& b) {\n";
        if (base) {
            cpp += "    if (!((const " + base->tname() + "&)a == (const " +
                   base->tname() + "&)b)) return false;\n";
        }
        for (auto var : vars) {
            cpp += "    if(!(a." + var->name + " == b." + var->name +
                   ")) return false;\n";
        }
        cpp += "    return true;\n";
        cpp += "}\n\n";
        return cpp;
    }
    virtual std::string to_destructor_func() {
        auto cpp = std::string();
        auto has_destructor = false;
        cpp += "\n /// destructor\n";
        cpp += "~" + tname() + "() {\n";
        for (auto var : vars) {
            if (var->type->is_struct()) {
                cpp += "if(" + var->name + ") delete " + var->name + ";\n";
                has_destructor = true;
            }
            if (var->type->is_array()) {
                if (!((ArrayType*)(var->type))->items->is_struct()) continue;
                cpp += "for(auto e : " + var->name + ") if(e) delete e;\n";
                has_destructor = true;
            }
        }
        cpp += "}";
        return (has_destructor) ? cpp : "";
    }
    virtual std::string to_default(const json& js) { return "nullptr"; }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname;
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return "!" + varname;
    }
};

struct DictType : Type {
    Type* items = nullptr;
    DictType(Type* items) : items(items) {}

    virtual std::string to_invalid_check() { return ".empty()"; }
    virtual std::string tname() {
        return "std::map<std::string, " + items->rname() + ">";
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return "!" + varname + ".empty()";
    }
    virtual std::string to_invalid_value(
        const std::string& varname, const json& def) {
        return varname + ".empty()";
    }
};

template <typename T>
struct EnumType : Type {
    std::string name = "";
    Type* base = nullptr;
    Type* parent = nullptr;
    std::string help = "";
    std::vector<std::pair<std::string, T>> labels = {};

    EnumType(const std::string& name, Type* base, Type* parent)
        : name(name), base(base), parent(parent) {}

    virtual std::string tname() { return name; }
    virtual std::string fname() {
        if (parent) return parent->fname() + "::" + tname();
        return tname();
    }

    static std::string _to_type_label(int v, int idx) {
        return std::to_string(v);
    }
    static std::string _to_type_label(const std::string& s, int idx) {
        return std::to_string(idx);
    }

    static std::string _to_parse_func_label(int i) { return std::to_string(i); }
    static std::string _to_parse_func_label(const std::string& s) {
        return "\"" + s + "\"";
    }

    virtual std::string to_type() {
        auto cpp = std::string();
        cpp += "///\n/// " + help + "\n///\n";
        cpp += "enum struct " + tname() + " {\n";
        auto count = 0;
        cpp += "/// not set\n";
        cpp += "NotSet = -1,\n";
        for (auto label : labels) {
            cpp += "/// " + label.first + "\n";
            cpp += label.first + " = " + _to_type_label(label.second, count++) +
                   ",\n";
        }
        cpp += "};\n\n";
        return cpp;
    }
    virtual std::string to_default(const json& js) {
        if (js.empty()) { return name + "::" + "NotSet"; }
        auto val = js.get<T>();
        for (auto label : labels)
            if (label.second == val) return name + "::" + label.first;
        assert(false);
        return "";
    }
    virtual std::string to_defined_value(
        const std::string& varname, const json& def) {
        return varname + " != " + tname() + "::NotSet";
    }
    virtual std::string to_parse_func() {
        auto cpp = std::string();
        cpp += "// Parse a " + fname() + " enum\n";
        cpp += "static bool parse(" + fname() +
               "& val, const json& js, parse_stack& err) {\n";
        cpp += "    static std::map<" + base->tname() + ", " + fname() +
               "> table = {\n";
        for (auto label : labels) {
            cpp += "        { " + _to_parse_func_label(label.second) + ", " +
                   fname() + "::" + label.first + " },\n";
        }
        cpp += "    };\n";
        cpp += "    auto v = " + base->tname() + "();\n";
        cpp += "    parse(v, js, err);\n";
        cpp += "    if (table.find(v) == table.end()) return false;\n";
        cpp += "    val = table[v];\n";
        cpp += "    return true;\n";
        cpp += "}\n\n";
        return cpp;
    }
    virtual std::string to_dump_func() {
        auto cpp = std::string();
        cpp += "// Converts a " + fname() + " enum to JSON\n";
        cpp += "static void dump(const " + fname() +
               "& val, json& js, parse_stack& err) {\n";
        cpp += "    static std::map<" + fname() + ", " + base->tname() +
               "> table = {\n";
        for (auto label : labels) {
            cpp += "        { " + fname() + "::" + label.first + ", " +
                   _to_parse_func_label(label.second) + " },\n";
        }
        cpp += "    };\n";
        cpp += "    auto v = table.at(val);\n";
        cpp += "    dump(v, js, err);\n";
        cpp += "}\n\n";
        return cpp;
    }
};

using IntEnumType = EnumType<int>;
using StringEnumType = EnumType<std::string>;

// clang-format off
auto filenames = std::vector<std::string>{
    "glTFProperty.schema.json",
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
};
// clang-format on

Type* schema2type(const std::vector<Type*>& types, const std::string& name) {
    for (auto type : types)
        if (type->sname() == name) return type;
    return nullptr;
}

StructType* schema2stype(
    const std::vector<Type*>& types, const std::string& name) {
    return (StructType*)schema2type(types, name);
}

json load_json(const std::string& filename) {
    auto stream = std::ifstream(filename);
    if (!stream) {
        printf("cannot open schema %s\n", filename.c_str());
        exit(1);
    }
    json js;
    stream >> js;
    return js;
}

Type* varschema2type(std::vector<Type*>& types, const json& schema,
    const std::string& varname, StructType* parent) {
    if (schema.count("type")) {
        auto& type = schema.at("type");
        if (type == "array") {
            auto nitems = 0;
            if (schema.count("minItems") && schema.count("maxItems")) {
                if (schema["minItems"] == schema["maxItems"])
                    nitems = schema["minItems"];
            }
            auto itype =
                varschema2type(types, schema.at("items"), varname, parent);
            if (itype->is_id())
                itype = new TypedIdType(
                    schemaname2typename(schema.at("gltf_id_type")));
            if (!nitems || itype->tname() != "float")
                return new ArrayType(itype);
            if (varname == "rotation")
                return new FixedArrayType("ym::quat4f", 4, 0);
            if (varname == "matrix")
                return new FixedArrayType("ym::mat4f", 4, 4);
            if (nitems == 2) return new FixedArrayType("ym::vec2f", 2, 0);
            if (nitems == 3) return new FixedArrayType("ym::vec3f", 3, 0);
            if (nitems == 4) return new FixedArrayType("ym::vec4f", 4, 0);
            assert(false);
            return nullptr;
        } else if (type == "object") {
            auto sname = filename2schemaname(
                schema.at("additionalProperties").at("$ref"));
            auto itype = schema2type(types, sname);
            if (itype->is_id())
                itype = new TypedIdType(
                    schemaname2typename(schema.at("gltf_id_type")));
            auto atype = new DictType(itype);
            return atype;
        } else {
            return schema2type(types, type);
        }
    } else if (schema.count("allOf")) {
        auto sname = filename2schemaname(schema["allOf"][0]["$ref"]);
        auto type = schema2type(types, sname);
        if (type->is_id()) {
            return new TypedIdType(
                schemaname2typename(schema.at("gltf_id_type")));
        } else {
            return type;
        }
    } else if (schema.count("$ref")) {
        auto sname = filename2schemaname(schema["$ref"]);
        return schema2type(types, sname);
    } else if (schema.count("anyOf")) {
        auto elems = schema.at("anyOf");
        auto ebase = elems.back().at("type").get<std::string>();
        if (ebase == "string") {
            auto etype =
                new StringEnumType(parent->tname() + varname2enumname(varname),
                    schema2type(types, "string"), nullptr);
            etype->help =
                "Values for " + parent->tname() + "::" + varname + ".";
            for (auto e : elems) {
                if (!e.count("enum")) continue;
                auto value = e.at("enum").at(0).get<std::string>();
                auto label = string2enumvalue(value);
                etype->labels += {label, value};
            }
            types += etype;
            return etype;
        } else {
            auto etype =
                new IntEnumType(parent->tname() + varname2enumname(varname),
                    schema2type(types, "integer"), nullptr);
            etype->help =
                "Enum values for " + parent->tname() + "::" + varname + ".";
            for (auto e : elems) {
                if (!e.count("enum")) continue;
                auto label = string2enumvalue(e.at("description"));
                auto value = e.at("enum").at(0).get<int>();
                etype->labels += {label, value};
            }

            types += etype;
            return etype;
        }
    }
    throw std::runtime_error("missing type");
    return nullptr;
}

void schemas2types(
    std::vector<Type*>& types, const std::string& schemadir, bool subtypes) {
    for (auto&& filename : filenames) {
        printf("parsing %s\n", filename.c_str());
        auto sname = filename2schemaname(filename);
        auto tname = schemaname2typename(sname);
        auto schema = load_json(schemadir + filename);
        auto type = new StructType(sname, tname);
        if (schema.count("description"))
            type->help = schema["description"].get<std::string>();
        if (schema.count("allOf")) {
            auto basename = filename2schemaname(schema["allOf"][0]["$ref"]);
            auto base = (StructType*)schema2type(types, basename);
            if (base) {
                type->base = base;
                base->subtypes += type;
            }
        }
        auto props = schema["properties"];
        auto required =
            (schema.count("required")) ? schema["required"] : json();
        for (auto var_it = props.begin(); var_it != props.end(); ++var_it) {
            auto var_schema = var_it.value();
            if (var_it.value().empty()) continue;
            auto req = false;
            for (auto req_it = required.begin(); req_it != required.end();
                 ++req_it) {
                if (var_it.key() == *req_it) req = true;
            }
            auto var = new Variable(var_it.key(),
                varschema2type(types, var_schema, var_it.key(), type));
            var->required = req;
            if (var_schema.count("gltf_extension"))
                var->ext_name = var_schema["gltf_extension"].get<std::string>();
            if (var_schema.count("description"))
                var->help = var_schema["description"].get<std::string>();
            if (var_schema.count("gltf_detailedDescription"))
                var->long_help =
                    var_schema["gltf_detailedDescription"].get<std::string>();
            if (var_schema.count("default"))
                var->defvalue = var_schema.at("default");
            type->vars += var;
        }
        types += type;
    }
}

std::string substitute(
    const std::string& txt, const std::string& code, const std::string& label) {
    auto lines = splitlines(txt, true);
    auto ret = std::string();
    auto skipping = false;
    for (auto& line : lines) {
        if (contains(line, "#codegen") && contains(line, "begin") &&
            contains(line, label)) {
            ret += line;
            ret += code;
            skipping = true;
        } else if (contains(line, "#codegen") && contains(line, "end") &&
                   contains(line, label)) {
            ret += line;
            skipping = false;
        } else {
            if (!skipping) ret += line;
        }
    }
    if (skipping) {
        printf("error parsing file: unmatched begin/end\n\n");
        exit(1);
    }
    return ret;
}

std::vector<Type*> init_types() {
    auto types = std::vector<Type*>{};
    types += new IntType();
    types += new FloatType();
    types += new BoolType();
    types += new StringType();
    types += new JsonType();
    types += new DictType(types.back());
    types += new IdType();
    types += new AliasType("extension", "extension_t", types[5]);
    types += new AliasType("extras", "extras_t", types[4]);
    return types;
}

void add_extra_vars(std::vector<Type*> types) {
    types += new ExternType("image.data", "image_data");
    auto image = (StructType*)schema2type(types, "image");
    image->extra_vars += new Variable("data", types.back());
    image->extra_vars.back()->help = "Image data if loaded.";
    types += new ExternType("buffer.data", "buffer_data");
    auto buffer = (StructType*)schema2type(types, "buffer");
    buffer->extra_vars += new Variable("data", types.back());
    buffer->extra_vars.back()->help = "Buffer data if loaded.";
}

int main(int argc, char** argv) {
    auto parser =
        make_parser(argc, argv, "generates gltf c++ code from json schema");
    auto schemadir = parse_opt<std::string>(parser, "--schemadir", "",
        "schema directory", "tools/gltf/schema/", false);
    auto type_filename = parse_opt<std::string>(parser, "--type-filename", "",
        "type filename", "yocto/yocto_gltf.h", false);
    auto func_filename = parse_opt<std::string>(parser, "--func-filename", "",
        "func filename", "yocto/yocto_gltf.cpp", false);
    auto standalone = parse_flag(parser, "--standalone", "", "standalone code");
    auto no_format =
        parse_flag(parser, "--no-format", "", "do not run clang format");
    auto subtypes = parse_flag(parser, "--subtypes", "-s", "use subtypes");
    check_parser(parser);

    auto types = init_types();
    schemas2types(types, schemadir, subtypes);
    add_extra_vars(types);

    // generate type code
    auto type_cpp = std::string();
    type_cpp += "// forward decalration\n";
    for (auto type : types)
        if (type->is_struct()) type_cpp += "struct " + type->tname() + ";\n";
    type_cpp += "\n";
    for (auto type : types) type_cpp += type->to_type();
    if (!standalone) {
        auto code = load_txtfile(type_filename);
        type_cpp = substitute(code, type_cpp, "type");
    }
    save_txtfile(type_filename, type_cpp);

    // generate func code
    auto func_cpp = std::string();
    func_cpp += "\n" + Type::to_parse_support_func() + "\n";
    for (auto type : types) func_cpp += type->to_parse_func();
    // for (auto type : types) func_cpp += type->to_equal_func();
    func_cpp += "\n" + Type::to_dump_support_func() + "\n";
    for (auto type : types) func_cpp += type->to_dump_func();
    func_cpp += "\n" + Type::to_validate_support_func() + "\n";
    for (auto type : types) func_cpp += type->to_validate_func();
    if (!standalone) {
        auto code = load_txtfile(func_filename);
        func_cpp = substitute(code, func_cpp, "func");
    }
    save_txtfile(func_filename, func_cpp);

    // calls clang format
    if (!no_format) {
        system(("/usr/local/bin/clang-format -i -style=file " + type_filename)
                   .c_str());
        system(("/usr/local/bin/clang-format -i -style=file " + func_filename)
                   .c_str());
    }
}
