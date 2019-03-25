//
// Implementation for Yocto/Pbrt loader.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_pbrt.h"

#include <memory>
#include <string_view>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF LOW LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Token stream
struct pbrt_token_stream {
    string      buffer;
    string_view str;
};

// Param type
enum struct pbrt_param_type { int_t, float_t, bool_t, string_t, texture_t, rgb_t };

static inline bool is_alpha(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}
static inline bool is_digit(char c) { return c >= '0' && c <= '9'; }

// skip white space or comment
void skip_whitespace_or_comment(pbrt_token_stream& stream) {
    auto& str = stream.str;
    if (str.empty()) return;
    while (!str.empty() && (std::isspace(str.front()) || str.front() == '#')) {
        if (str.front() == '#') {
            auto pos = str.find('\n');
            if (pos != string_view::npos) {
                str.remove_prefix(pos);
            } else {
                str.remove_prefix(str.length());
            }
        } else {
            auto pos = str.find_first_not_of(" \t\n\r");
            if (pos == string_view::npos) {
                str.remove_prefix(str.length());
            } else {
                str.remove_prefix(pos);
            }
        }
    }
}

// parse a quoted string
static inline void parse_value(pbrt_token_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.front() != '"') {
        throw pbrtio_error("bad string");
    }
    str.remove_prefix(1);
    auto pos = str.find('"');
    if (pos == string_view::npos) {
        throw pbrtio_error("bad string");
    }
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
}

// parse a quoted string
static inline void parse_command(pbrt_token_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (!std::isalpha((int)str.front())) {
        throw pbrtio_error("bad command");
    }
    auto pos = str.find_first_not_of(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    if (pos == string_view::npos) {
        value.assign(str);
        str.remove_prefix(str.size());
    } else {
        value.assign(str.substr(0, pos));
        str.remove_prefix(pos + 1);
    }
}

// parse a number
static inline void parse_value(pbrt_token_stream& stream, float& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw pbrtio_error("number expected");
    auto next = (char*)nullptr;
    value     = strtof(str.data(), &next);
    if (str.data() == next) throw pbrtio_error("number expected");
    str.remove_prefix(next - str.data());
}

// parse a number
static inline void parse_value(pbrt_token_stream& stream, int& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw pbrtio_error("number expected");
    auto next = (char*)nullptr;
    value     = strtol(str.data(), &next, 10);
    if (str.data() == next) throw pbrtio_error("number expected");
    str.remove_prefix(next - str.data());
}
static inline void parse_value(pbrt_token_stream& stream, bool& value) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    if (value_name == "true") {
        value = true;
    } else if (value_name == "false") {
        value = false;
    } else {
        throw pbrtio_error("expected boolean");
    }
}

// parse a vec type
template <typename T, int N>
static inline void parse_value(pbrt_token_stream& stream, vec<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}
template <typename T, int N, int M>
static inline void parse_value(pbrt_token_stream& stream, mat<T, N, M>& value) {
    for (auto i = 0; i < M; i++) parse_value(stream, value[i]);
}
template <typename T>
static inline void parse_value(pbrt_token_stream& stream, bbox<T, 2>& value) {
    parse_value(stream, value[0][0]);
    parse_value(stream, value[1][0]);
    parse_value(stream, value[0][1]);
    parse_value(stream, value[1][1]);
}

// Check if empty
inline bool is_empty(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return stream.str.empty();
}

// Check if empty
inline bool is_string(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '"';
}
inline bool is_open_bracket(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '[';
}
inline bool is_close_bracket(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == ']';
}

// parse a quoted string
static inline void parse_nametype(
    pbrt_token_stream& stream, string& name, pbrt_param_type& type) {
    static const auto types = unordered_map<string, pbrt_param_type>{
        {"integer", pbrt_param_type::int_t},
        {"float", pbrt_param_type::float_t},
        {"boolean", pbrt_param_type::bool_t},
        {"string", pbrt_param_type::string_t},
        {"texture", pbrt_param_type::texture_t},
        {"rgb", pbrt_param_type::rgb_t},
        {"color", pbrt_param_type::rgb_t},
    };
    auto value = ""s;
    parse_value(stream, value);
    auto str  = string_view{value};
    auto pos1 = str.find(' ');
    if (pos1 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    auto type_name = string(str.substr(0, pos1));
    str.remove_prefix(pos1);
    auto pos2 = str.find_first_not_of(' ');
    if (pos2 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    str.remove_prefix(pos2);
    name = string(str);
    try {
        type = types.at(type_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("unknown type "s + type_name);
    }
}

static inline void skip_open_bracket(pbrt_token_stream& stream) {
    if (!is_open_bracket(stream)) throw pbrtio_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}
static inline void skip_close_bracket(pbrt_token_stream& stream) {
    if (!is_close_bracket(stream)) throw pbrtio_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}

static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, int& value) {
    if (type != pbrt_param_type::int_t) throw pbrtio_error("expected int");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, float& value) {
    if (type != pbrt_param_type::float_t) throw pbrtio_error("expected int");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, bbox2i& value) {
    if (type != pbrt_param_type::int_t) throw pbrtio_error("expected bbox2i");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(pbrt_token_stream& stream, pbrt_param_type type, bbox2f& value) {
    if (type != pbrt_param_type::float_t) throw pbrtio_error("expected bbox2f");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, vec3f& value) {
    if (type != pbrt_param_type::rgb_t) throw pbrtio_error("expected vec4i");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, mat4f& value) {
    if (type != pbrt_param_type::float_t) throw pbrtio_error("expected vec4i");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, bool& value) {
    if (type != pbrt_param_type::bool_t) throw pbrtio_error("expected int");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, string& value) {
    if (type != pbrt_param_type::string_t && type != pbrt_param_type::texture_t) throw pbrtio_error("expected string");
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}

template<typename T>
static inline void parse_param(
    pbrt_token_stream& stream, pbrt_param_type type, pbrt_textured<T>& value) {
    if (type == pbrt_param_type::texture_t) {
        parse_param(stream, type, value.texture);
    } else {
        parse_param(stream, type, value.value);
    }
}

// Load a token stream
static void load_token_stream(
    const string& filename, pbrt_token_stream& stream) {
    load_text(filename, stream.buffer);
    stream.str = stream.buffer;
}

// operations on token stacks
static void init_token_streams(vector<pbrt_token_stream>& streams) {
    streams.reserve(100);
}
static void load_token_stream(
    const string& filename, vector<pbrt_token_stream>& streams) {
    streams.emplace_back();
    load_token_stream(filename, streams.back());
}

// Skip whitespace
inline void skip_whitespace_or_comment(vector<pbrt_token_stream>& streams) {
    if (streams.empty()) return;
    while (!streams.empty()) {
        skip_whitespace_or_comment(streams.back());
        if (is_empty(streams.back())) {
            streams.pop_back();
        } else {
            break;
        }
    }
}

// Check if empty
inline bool is_empty(vector<pbrt_token_stream>& streams) {
    skip_whitespace_or_comment(streams);
    return streams.empty();
}
inline bool is_param(vector<pbrt_token_stream>& streams) {
    skip_whitespace_or_comment(streams);
    return !streams.empty() && is_string(streams.back());
}

// Parse values
static inline void parse_command(
    vector<pbrt_token_stream>& streams, string& value) {
    return parse_command(streams.back(), value);
}
template <typename T>
static inline void parse_value(vector<pbrt_token_stream>& streams, T& value) {
    return parse_value(streams.back(), value);
}
static inline void parse_nametype(
    vector<pbrt_token_stream>& streams, string& name, pbrt_param_type& type) {
    return parse_nametype(streams.back(), name, type);
}
template <typename T>
static inline void parse_param(
    vector<pbrt_token_stream>& streams, pbrt_param_type ptype, T& value) {
    return parse_param(streams.back(), ptype, value);
}
template <typename T>
static inline void parse_optional_param(vector<pbrt_token_stream>& streams,
    const string& pname, pbrt_param_type ptype, const char* name, T& value) {
    if (pname != name) return;
    return parse_param(streams, ptype, value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Parse Integrator
void parse_pbrt_integrator(
    vector<pbrt_token_stream>& streams, pbrt_integrator& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    if (type == "path") {
        auto tvalue = pbrt_integrator_path{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "rrthreshold") {
                parse_param(streams, ptype, tvalue.rrthreshold);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
            // parse_optional_param(streams, "lightsamplestrategy",
            // tvalue.lightsamplestrategy); // TODO: enums
        }
        value = tvalue;
    } else if (type == "directlighting") {
        auto tvalue = pbrt_integrator_directlighting{};
        throw pbrtio_error("unsupported Integrator " + type);
        value = tvalue;
    } else if (type == "bdpt") {
        auto tvalue = pbrt_integrator_bdpt{};
        throw pbrtio_error("unsupported Integrator " + type);
        value = tvalue;
    } else if (type == "mlt") {
        auto tvalue = pbrt_integrator_mlt{};
        throw pbrtio_error("unsupported Integrator " + type);
        value = tvalue;
    } else if (type == "sppm") {
        auto tvalue = pbrt_integrator_sppm{};
        throw pbrtio_error("unsupported Integrator " + type);
        value = tvalue;
    } else if (type == "whitted") {
        auto tvalue = pbrt_integrator_whitted{};
        throw pbrtio_error("unsupported Integrator " + type);
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Integrator " + type);
    }
}

// Parse Sampler
void parse_pbrt_sampler(
    vector<pbrt_token_stream>& streams, pbrt_sampler& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    if (type == "random") {
        auto tvalue = pbrt_sampler_random{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(streams, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "halton") {
        auto tvalue = pbrt_sampler_halton{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(streams, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sobol") {
        auto tvalue = pbrt_sampler_sobol{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(streams, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "zerotwosequence") {
        auto tvalue = pbrt_sampler_zerotwosequence{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(streams, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "maxmindist") {
        auto tvalue = pbrt_sampler_maxmindist{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(streams, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "stratified") {
        auto tvalue = pbrt_sampler_stratified{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xsamples") {
                parse_param(streams, ptype, tvalue.xsamples);
            } else if (pname == "ysamples") {
                parse_param(streams, ptype, tvalue.ysamples);
            } else if (pname == "jitter") {
                parse_param(streams, ptype, tvalue.jitter);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Sampler " + type);
    }
}

// Parse Filter
void parse_pbrt_filter(vector<pbrt_token_stream>& streams, pbrt_filter& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    if (type == "box") {
        auto tvalue = pbrt_filter_box{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xwidth") {
                parse_param(streams, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(streams, ptype, tvalue.ywidth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "gaussian") {
        auto tvalue = pbrt_filter_gaussian{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xwidth") {
                parse_param(streams, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(streams, ptype, tvalue.ywidth);
            } else if (pname == "alpha") {
                parse_param(streams, ptype, tvalue.alpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mitchell") {
        auto tvalue = pbrt_filter_mitchell{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xwidth") {
                parse_param(streams, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(streams, ptype, tvalue.ywidth);
            } else if (pname == "B") {
                parse_param(streams, ptype, tvalue.B);
            } else if (pname == "C") {
                parse_param(streams, ptype, tvalue.C);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sinc") {
        auto tvalue = pbrt_filter_sinc{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xwidth") {
                parse_param(streams, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(streams, ptype, tvalue.ywidth);
            } else if (pname == "tau") {
                parse_param(streams, ptype, tvalue.tau);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "triangle") {
        auto tvalue = pbrt_filter_triangle{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xwidth") {
                parse_param(streams, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(streams, ptype, tvalue.ywidth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown PixelFilter " + type);
    }
}

// Parse Filter
void parse_pbrt_camera(vector<pbrt_token_stream>& streams, pbrt_camera& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    if (type == "perspective") {
        auto tvalue = pbrt_camera_perspective{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "fov") {
                parse_param(streams, ptype, tvalue.fov);
            } else if (pname == "frameaspectratio") {
                parse_param(streams, ptype, tvalue.frameaspectratio);
            } else if (pname == "lensradius") {
                parse_param(streams, ptype, tvalue.lensradius);
            } else if (pname == "focaldistance") {
                parse_param(streams, ptype, tvalue.focaldistance);
            } else if (pname == "screenwindow") {
                parse_param(streams, ptype, tvalue.screenwindow);
            } else if (pname == "shutteropen") {
                parse_param(streams, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(streams, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "orthographic") {
        auto tvalue = pbrt_camera_orthographic{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "frameaspectratio") {
                parse_param(streams, ptype, tvalue.frameaspectratio);
            } else if (pname == "lensradius") {
                parse_param(streams, ptype, tvalue.lensradius);
            } else if (pname == "focaldistance") {
                parse_param(streams, ptype, tvalue.focaldistance);
            } else if (pname == "screenwindow") {
                parse_param(streams, ptype, tvalue.screenwindow);
            } else if (pname == "shutteropen") {
                parse_param(streams, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(streams, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "environment") {
        auto tvalue = pbrt_camera_environment{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "shutteropen") {
                parse_param(streams, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(streams, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "realistic") {
        auto tvalue = pbrt_camera_realistic{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "lensfile") {
                parse_param(streams, ptype, tvalue.lensfile);
            } else if (pname == "aperturediameter") {
                parse_param(streams, ptype, tvalue.aperturediameter);
            } else if (pname == "focusdistance") {
                parse_param(streams, ptype, tvalue.focusdistance);
            } else if (pname == "simpleweighting") {
                parse_param(streams, ptype, tvalue.simpleweighting);
            } else if (pname == "shutteropen") {
                parse_param(streams, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(streams, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Filter
void parse_pbrt_material(
    vector<pbrt_token_stream>& streams, pbrt_material& value, bool named) {
    static int material_id = 0;
    auto name = ""s, type = ""s;
    if(!named) {
        parse_value(streams, type);
        name = "unnamed_material_" + std::to_string(material_id++);
    } else {
        parse_value(streams, name);
        auto ptname = ""s;
        auto pttype = pbrt_param_type::string_t;
        parse_nametype(streams, ptname, pttype);
        if(ptname == "type") {
            parse_param(streams, pttype, type);
        } else {
            throw pbrtio_error("expected material type");
        }
    }
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    struct pbrt_material_plastic {
        string               name           = "";
        pbrt_textured<vec3f> Kd             = {0.25f, 0.25f, 0.25f};
        pbrt_textured<vec3f> Ks             = {0.25f, 0.25f, 0.25f};
        pbrt_textured<float> roughness      = 0.1f;
        bool                 remaproughness = true;
    };
    struct pbrt_material_metal {
        string               name = "";
        pbrt_textured<vec3f> eta  = {
            0.2004376970f, 0.9240334304f, 1.1022119527f};
        pbrt_textured<vec3f> k = {3.9129485033f, 2.4528477015f, 2.1421879552f};
        pbrt_textured<float> roughness      = 0.01;
        pbrt_textured<float> uroughness     = 0.01;
        pbrt_textured<float> vroughness     = 0.01;
        bool                 remaproughness = true;
    };
    struct pbrt_material_glass {
        string               name           = "";
        pbrt_textured<vec3f> Kr             = {1, 1, 1};
        pbrt_textured<vec3f> Kt             = {1, 1, 1};
        pbrt_textured<float> eta            = 1;
        pbrt_textured<float> uroughness     = 0;
        pbrt_textured<float> vroughness     = 0;
        bool                 remaproughness = true;
    };
    struct pbrt_material_translucent {
        string               name           = "";
        pbrt_textured<vec3f> Kd             = {0, 0, 0};
        pbrt_textured<vec3f> Ks             = {0, 0, 0};
        pbrt_textured<vec3f> reflect        = {0, 0, 0};
        pbrt_textured<vec3f> transmit       = {0, 0, 0};
        pbrt_textured<float> roughness      = 0;
        bool                 remaproughness = true;
    };
    struct pbrt_material_uber {
        string               name           = "";
        pbrt_textured<vec3f> Kd             = {0, 0, 0};
        pbrt_textured<vec3f> Ks             = {0, 0, 0};
        pbrt_textured<vec3f> Kr             = {0, 0, 0};
        pbrt_textured<vec3f> Kt             = {0, 0, 0};
        pbrt_textured<float> roughness      = 0;
        pbrt_textured<float> uroughness     = 0;
        pbrt_textured<float> vroughness     = 0;
        pbrt_textured<float> eta            = 1;
        pbrt_textured<vec3f> opacity        = {1, 1, 1};
        bool                 remaproughness = true;
    };
    struct pbrt_material_disney {
        string               name            = "";
        pbrt_textured<vec3f> color           = {0.5f, 0.5f, 0.5f};
        pbrt_textured<float> anisotropic     = 0;
        pbrt_textured<float> clearcoat       = 0;
        pbrt_textured<float> clearcoatgloss  = 1;
        pbrt_textured<float> eta             = 1.5f;
        pbrt_textured<float> metallic        = 0;
        pbrt_textured<float> roughness       = 0.5f;
        pbrt_textured<vec3f> scatterdistance = {0, 0, 0};
        pbrt_textured<float> sheen           = 0;
        pbrt_textured<float> sheentint       = 0.5;
        pbrt_textured<float> spectrans       = 0;
        pbrt_textured<float> speculartint    = 0;
        bool                 thin            = false;
        pbrt_textured<vec3f> difftrans       = {1, 1, 1};
        pbrt_textured<vec3f> flatness        = {0, 0, 0};
    };
    struct pbrt_material_fourier {
        string name     = "";
        string bsdffile = "";
    };
    struct pbrt_material_hair {
        string               name        = "";
        pbrt_textured<vec3f> sigma_a     = {0, 0, 0};  // TODO: missing default
        pbrt_textured<vec3f> color       = {0, 0, 0};  // TODO: missing default
        pbrt_textured<float> eumelanin   = 0;          // TODO: missing default
        pbrt_textured<float> pheomelanin = 0;          // TODO: missing default
        pbrt_textured<float> eta         = 1.55f;
        pbrt_textured<float> beta_m      = 0.3f;
        pbrt_textured<float> beta_n      = 0.3f;
        pbrt_textured<float> alpha       = 2;
    };
    struct pbrt_material_kdsubsurface {
        string               name           = "";
        pbrt_textured<vec3f> Kd             = {0, 0, 0};
        pbrt_textured<float> mfp            = 1;
        pbrt_textured<float> eta            = 1;
        pbrt_textured<vec3f> Kr             = {1, 1, 1};
        pbrt_textured<vec3f> Kt             = {1, 1, 1};
        pbrt_textured<float> uroughness     = 0;
        pbrt_textured<float> vroughness     = 0;
        bool                 remaproughness = true;
    };
    struct pbrt_material_mix {
        string               name           = "";
        pbrt_textured<vec3f> amount         = {0, 0, 0};
        string               namedmaterial1 = "";
        string               namedmaterial2 = "";
    };
    struct pbrt_material_substrate {
        string               name           = "";
        pbrt_textured<vec3f> Kd             = {0, 0, 0};
        pbrt_textured<vec3f> Ks             = {0, 0, 0};
        pbrt_textured<float> uroughness     = 0;
        pbrt_textured<float> vroughness     = 0;
        bool                 remaproughness = true;
    };
    struct pbrt_material_subsurface {
        string               name           = "";
        string               name_          = "";
        pbrt_textured<vec3f> sigma_a        = {.0011, .0024, .014};
        pbrt_textured<vec3f> sigma_prime_s  = {2.55, 3.12, 3.77};
        float                scale          = 1;
        pbrt_textured<float> eta            = 1;
        pbrt_textured<vec3f> Kr             = {1, 1, 1};
        pbrt_textured<vec3f> Kt             = {1, 1, 1};
        pbrt_textured<float> uroughness     = 0;
        pbrt_textured<float> vroughness     = 0;
        bool                 remaproughness = true;
    };
    if (type == "matte") {
        auto tvalue = pbrt_material_matte{};
        tvalue.name = name;
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "sigma") {
                parse_param(streams, ptype, tvalue.sigma);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mirror") {
        auto tvalue = pbrt_material_mirror{};
        tvalue.name = name;
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Filter
void parse_pbrt_film(vector<pbrt_token_stream>& streams, pbrt_film& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s;
    auto ptype = pbrt_param_type::int_t;
    if (type == "image") {
        auto tvalue = pbrt_film_image{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "xresolution") {
                parse_param(streams, ptype, tvalue.xresolution);
            } else if (pname == "yresolution") {
                parse_param(streams, ptype, tvalue.yresolution);
            } else if (pname == "yresolution") {
                parse_param(streams, ptype, tvalue.yresolution);
            } else if (pname == "cropwindow") {
                parse_param(streams, ptype, tvalue.cropwindow);
            } else if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "maxsampleluminance") {
                parse_param(streams, ptype, tvalue.maxsampleluminance);
            } else if (pname == "diagonal") {
                parse_param(streams, ptype, tvalue.diagonal);
            } else if (pname == "filename") {
                parse_param(streams, ptype, tvalue.filename);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Load pbrt scene
void load_pbrt(const string& filename, const pbrt_callbacks& cb,
    const load_pbrt_options& options) {
    // start laoding files
    auto streams = vector<pbrt_token_stream>{};
    init_token_streams(streams);
    load_token_stream(filename, streams);

    // parse command by command
    auto cmd = ""s;
    while (!is_empty(streams)) {
        // get command
        parse_command(streams, cmd);
        if (cmd == "WorldBegin") {
        } else if (cmd == "WorldEnd") {
        } else if (cmd == "Transform") {
            auto xf = identity_mat4f;
            parse_param(streams, pbrt_param_type::float_t, xf);
        } else if (cmd == "Integrator") {
            auto value = pbrt_integrator{};
            parse_pbrt_integrator(streams, value);
        } else if (cmd == "Sampler") {
            auto value = pbrt_sampler{};
            parse_pbrt_sampler(streams, value);
        } else if (cmd == "PixelFilter") {
            auto value = pbrt_filter{};
            parse_pbrt_filter(streams, value);
        } else if (cmd == "Film") {
            auto value = pbrt_film{};
            parse_pbrt_film(streams, value);
        } else if (cmd == "Camera") {
            auto value = pbrt_camera{};
            parse_pbrt_camera(streams, value);
        } else if (cmd == "MakeNamedMaterial") {
            auto value = pbrt_material{};
            parse_pbrt_material(streams, value, true);
        } else if (cmd == "NamedMaterial") {
            auto name = ""s;
            parse_value(streams, name);
        } else {
            throw pbrtio_error("unknown command " + cmd);
        }
    }
}

}  // namespace yocto
