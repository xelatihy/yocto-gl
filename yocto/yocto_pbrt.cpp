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
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_bilerp::mapping_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_bilerp::mapping_type>{
            {"uv", pbrt_texture_bilerp::mapping_type::uv},
            {"spherical", pbrt_texture_bilerp::mapping_type::spherical},
            {"cylindrical", pbrt_texture_bilerp::mapping_type::cylindrical},
            {"planar", pbrt_texture_bilerp::mapping_type::planar},
        };
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("expected mapping_type");
    }
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_checkerboard::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_dots::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_imagemap::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_uv::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}

static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_checkerboard::aamode_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_checkerboard::aamode_type>{
            {"closedform", pbrt_texture_checkerboard::aamode_type::closedform},
            {"none", pbrt_texture_checkerboard::aamode_type::none},
        };
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("expected aamode_type");
    }
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_imagemap::wrap_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_imagemap::wrap_type>{
            {"repeat", pbrt_texture_imagemap::wrap_type::repeat},
            {"clamp", pbrt_texture_imagemap::wrap_type::clamp},
            {"black", pbrt_texture_imagemap::wrap_type::black},
        };
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("expected aamode_type");
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
    pbrt_token_stream& stream, string& name, string& type) {
    auto value = ""s;
    parse_value(stream, value);
    auto str  = string_view{value};
    auto pos1 = str.find(' ');
    if (pos1 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    type = string(str.substr(0, pos1));
    str.remove_prefix(pos1);
    auto pos2 = str.find_first_not_of(' ');
    if (pos2 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    str.remove_prefix(pos2);
    name = string(str);
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

template <typename T>
static inline void parse_param(pbrt_token_stream& stream, T& value) {
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}

template <typename T>
static inline void parse_param(pbrt_token_stream& stream, vector<T>& values) {
    skip_open_bracket(stream);
    values.clear();
    while (!is_close_bracket(stream)) {
        values.push_back({});
        parse_value(stream, values.back());
    }
    skip_close_bracket(stream);
}

template <typename T>
static inline bool is_type_compatible(const string& type) {
    if constexpr (std::is_same<T, int>::value) {
        return type == "integer";
    } else if constexpr (std::is_same<T, float>::value) {
        return type == "float";
    } else if constexpr (std::is_same<T, bool>::value) {
        return type == "boolean";
    } else if constexpr (std::is_same<T, string>::value) {
        return type == "string";
    } else if constexpr (std::is_same<T, vec2f>::value) {
        return type == "point2" || type == "vector2" || type == "float";
    } else if constexpr (std::is_same<T, vec3f>::value) {
        return type == "point3" || type == "vector3" || type == "normal3" ||
               type == "point" || type == "vector" || type == "normal" ||
               type == "float";
    } else if constexpr (std::is_same<T, vec3i>::value) {
        return type == "integer";
    } else if constexpr (std::is_same<T, bbox2i>::value) {
        return type == "integer";
    } else if constexpr (std::is_same<T, bbox2f>::value) {
        return type == "float";
    } else if constexpr (std::is_enum<T>::value) {
        return type == "string";
    } else {
        return false;
    }
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, T& value) {
    if (!is_type_compatible<T>(type)) {
        throw pbrtio_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, vector<T>& value) {
    if (!is_type_compatible<T>(type)) {
        throw pbrtio_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

static inline void parse_param(pbrt_token_stream& stream, const string& type,
    pbrt_textured<float>& value) {
    if (type == "texture") {
        parse_param(stream, type, value.texture);
    } else if (type == "float") {
        parse_param(stream, type, value.value);
    } else {
        throw pbrtio_error("incomparible textured type " + type);
    }
}

static inline void parse_param(pbrt_token_stream& stream, const string& type,
    pbrt_textured<vec3f>& value) {
    if (type == "texture") {
        parse_param(stream, value.texture);
    } else if (type == "rgb") {
        parse_param(stream, value.value);
    } else if (type == "spectrum") {
        printf("spectrum not supported well\n");
        value = {1,0,0};
        // throw pbrtio_error("spectrum not supported");
    } else {
        throw pbrtio_error("incomparible textured type " + type);
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
    vector<pbrt_token_stream>& streams, string& name, string& type) {
    return parse_nametype(streams.back(), name, type);
}
template <typename T>
static inline void parse_param(
    vector<pbrt_token_stream>& streams, const string& ptype, T& value) {
    return parse_param(streams.back(), ptype, value);
}
template <typename T>
static inline void parse_param(vector<pbrt_token_stream>& streams, T& value) {
    return parse_param(streams.back(), value);
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
    auto pname = ""s, ptype = ""s;
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
    auto pname = ""s, ptype = ""s;
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
    auto pname = ""s, ptype = ""s;
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
void parse_pbrt_film(vector<pbrt_token_stream>& streams, pbrt_film& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
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

// Parse Camera
void parse_pbrt_camera(vector<pbrt_token_stream>& streams, pbrt_camera& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
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

// Parse Texture
void parse_pbrt_texture(
    vector<pbrt_token_stream>& streams, pbrt_texture& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
    if (type == "constant") {
        auto tvalue = pbrt_texture_constant{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "value") {
                parse_param(streams, ptype, tvalue.value);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bilerp") {
        auto tvalue = pbrt_texture_bilerp{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "v00") {
                parse_param(streams, ptype, tvalue.v00);
            } else if (pname == "v01") {
                parse_param(streams, ptype, tvalue.v01);
            } else if (pname == "v10") {
                parse_param(streams, ptype, tvalue.v10);
            } else if (pname == "v11") {
                parse_param(streams, ptype, tvalue.v11);
            } else if (pname == "mapping") {
                parse_param(streams, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(streams, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(streams, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(streams, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(streams, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(streams, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "checkerboard") {
        auto tvalue = pbrt_texture_checkerboard{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "dimension") {
                parse_param(streams, ptype, tvalue.dimension);
            } else if (pname == "tex1") {
                parse_param(streams, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(streams, ptype, tvalue.tex2);
            } else if (pname == "aamode") {
                parse_param(streams, ptype, tvalue.aamode);
            } else if (pname == "mapping") {
                parse_param(streams, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(streams, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(streams, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(streams, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(streams, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(streams, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "dots") {
        auto tvalue = pbrt_texture_dots{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "inside") {
                parse_param(streams, ptype, tvalue.inside);
            } else if (pname == "outside") {
                parse_param(streams, ptype, tvalue.outside);
            } else if (pname == "mapping") {
                parse_param(streams, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(streams, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(streams, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(streams, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(streams, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(streams, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "imagemap") {
        auto tvalue = pbrt_texture_imagemap{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "filename") {
                parse_param(streams, ptype, tvalue.filename);
            } else if (pname == "wrap") {
                parse_param(streams, ptype, tvalue.wrap);
            } else if (pname == "maxanisotropy") {
                parse_param(streams, ptype, tvalue.maxanisotropy);
            } else if (pname == "trilinear") {
                parse_param(streams, ptype, tvalue.trilinear);
            } else if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "gamma") {
                parse_param(streams, ptype, tvalue.gamma);
            } else if (pname == "mapping") {
                parse_param(streams, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(streams, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(streams, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(streams, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(streams, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(streams, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_texture_mix{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "tex1") {
                parse_param(streams, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(streams, ptype, tvalue.tex2);
            } else if (pname == "amount") {
                parse_param(streams, ptype, tvalue.amount);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "scale") {
        auto tvalue = pbrt_texture_scale{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "tex1") {
                parse_param(streams, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(streams, ptype, tvalue.tex2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fbm") {
        auto tvalue = pbrt_texture_fbm{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "octaves") {
                parse_param(streams, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(streams, ptype, tvalue.roughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "wrinkled") {
        auto tvalue = pbrt_texture_wrinkled{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "octaves") {
                parse_param(streams, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(streams, ptype, tvalue.roughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "windy") {
        auto tvalue = pbrt_texture_windy{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "") {
                // TODO: missing params
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "marble") {
        auto tvalue = pbrt_texture_marble{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "octaves") {
                parse_param(streams, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(streams, ptype, tvalue.roughness);
            } else if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "variation") {
                parse_param(streams, ptype, tvalue.variation);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uv") {
        auto tvalue = pbrt_texture_uv{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "mapping") {
                parse_param(streams, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(streams, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(streams, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(streams, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(streams, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(streams, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Material " + type);
    }
}

// Parse Material
void parse_pbrt_material(
    vector<pbrt_token_stream>& streams, pbrt_material& value, bool named) {
    auto type = ""s;
    if (!named) {
        parse_value(streams, type);
    } else {
        auto ptname = ""s, pttype = ""s;
        parse_nametype(streams, ptname, pttype);
        if (ptname == "type") {
            parse_param(streams, pttype, type);
        } else {
            throw pbrtio_error("expected material type");
        }
    }
    auto pname = ""s, ptype = ""s;
    if (type == "matte") {
        auto tvalue = pbrt_material_matte{};
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
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "plastic") {
        auto tvalue = pbrt_material_plastic{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(streams, ptype, tvalue.Ks);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "metal") {
        auto tvalue = pbrt_material_metal{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "k") {
                parse_param(streams, ptype, tvalue.k);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "glass") {
        auto tvalue = pbrt_material_glass{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(streams, ptype, tvalue.Kt);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "translucent") {
        auto tvalue = pbrt_material_translucent{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(streams, ptype, tvalue.Ks);
            } else if (pname == "reflect") {
                parse_param(streams, ptype, tvalue.reflect);
            } else if (pname == "transmit") {
                parse_param(streams, ptype, tvalue.transmit);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uber") {
        auto tvalue = pbrt_material_uber{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(streams, ptype, tvalue.Ks);
            } else if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(streams, ptype, tvalue.Kt);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "opacity") {
                parse_param(streams, ptype, tvalue.opacity);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disney") {
        auto tvalue = pbrt_material_disney{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "color") {
                parse_param(streams, ptype, tvalue.color);
            } else if (pname == "anisotropic") {
                parse_param(streams, ptype, tvalue.anisotropic);
            } else if (pname == "clearcoat") {
                parse_param(streams, ptype, tvalue.clearcoat);
            } else if (pname == "clearcoatgloss") {
                parse_param(streams, ptype, tvalue.clearcoatgloss);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "metallic") {
                parse_param(streams, ptype, tvalue.metallic);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else if (pname == "scatterdistance") {
                parse_param(streams, ptype, tvalue.scatterdistance);
            } else if (pname == "sheen") {
                parse_param(streams, ptype, tvalue.sheen);
            } else if (pname == "sheentint") {
                parse_param(streams, ptype, tvalue.sheentint);
            } else if (pname == "spectrans") {
                parse_param(streams, ptype, tvalue.spectrans);
            } else if (pname == "thin") {
                parse_param(streams, ptype, tvalue.thin);
            } else if (pname == "difftrans") {
                parse_param(streams, ptype, tvalue.difftrans);
            } else if (pname == "flatness") {
                parse_param(streams, ptype, tvalue.flatness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hair") {
        auto tvalue = pbrt_material_hair{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "color") {
                parse_param(streams, ptype, tvalue.color);
            } else if (pname == "sigma_a") {
                parse_param(streams, ptype, tvalue.sigma_a);
            } else if (pname == "eumelanin") {
                parse_param(streams, ptype, tvalue.eumelanin);
            } else if (pname == "pheomelanin") {
                parse_param(streams, ptype, tvalue.pheomelanin);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "beta_m") {
                parse_param(streams, ptype, tvalue.beta_m);
            } else if (pname == "beta_n") {
                parse_param(streams, ptype, tvalue.beta_n);
            } else if (pname == "alpha") {
                parse_param(streams, ptype, tvalue.alpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "kdsubsurface") {
        auto tvalue = pbrt_material_kdsubsurface{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(streams, ptype, tvalue.Kt);
            } else if (pname == "mfp") {
                parse_param(streams, ptype, tvalue.mfp);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_material_mix{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "amount") {
                parse_param(streams, ptype, tvalue.amount);
            } else if (pname == "namedmaterial1") {
                parse_param(streams, ptype, tvalue.namedmaterial1);
            } else if (pname == "namedmaterial2") {
                parse_param(streams, ptype, tvalue.namedmaterial2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fourier") {
        auto tvalue = pbrt_material_fourier{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "bsdffile") {
                parse_param(streams, ptype, tvalue.bsdffile);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "substrate") {
        auto tvalue = pbrt_material_substrate{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "Kd") {
                parse_param(streams, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(streams, ptype, tvalue.Ks);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "subsurface") {
        auto tvalue = pbrt_material_subsurface{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "name") {
                parse_param(streams, ptype, tvalue.name);
            } else if (pname == "sigma_a") {
                parse_param(streams, ptype, tvalue.sigma_a);
            } else if (pname == "sigma_prime_s") {
                parse_param(streams, ptype, tvalue.sigma_prime_s);
            } else if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "eta") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "Kr") {
                parse_param(streams, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(streams, ptype, tvalue.Kt);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(streams, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(streams, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(streams, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(streams, ptype, tvalue.remaproughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Material " + type);
    }
}

// Parse Shape
void parse_pbrt_shape(vector<pbrt_token_stream>& streams, pbrt_shape& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
    struct pbrt_shape_plymesh {
        string filename = {};
        // texture alpha
        // texture shadowalpha
    };
    struct pbrt_shape_curve {
        enum struct type_t { flat, ribbon, cylinder };
        enum struct basis_t { bezier, bspline };
        vector<vec3f> P          = {};
        basis_t       basis      = basis_t::bezier;
        int           degree     = 3;
        type_t        type       = type_t::flat;
        vector<vec3f> N          = {};
        float         width      = 1;
        float         width0     = 1;
        float         width1     = 1;
        int           splitdepth = 3;
    };
    struct pbrt_shape_loopsubdiv {
        int           levels  = 3;
        vector<int>   indices = {};
        vector<vec3f> P       = {};
    };
    struct pbrt_shape_nurbs {
        int           nu     = -1;
        int           nv     = -1;
        vector<float> uknots = {};
        vector<float> vknots = {};
        float         u0     = -1;
        float         v0     = -1;
        float         u1     = -1;
        float         v1     = -1;
        vector<vec3f> P      = {};
        vector<float> Pw     = {};
    };
    struct pbrt_shape_sphere {
        float radius = 1;
        float zmin   = -radius;
        float zmax   = radius;
        float phimax = 360;
    };
    struct pbrt_shape_disk {
        float height      = 0;
        float radius      = 1;
        float innerradius = 0;
        float phimax      = 360;
    };
    struct pbrt_shape_cone {
        float radius = 1;
        float height = 1;
        float phimax = 360;
    };
    struct pbrt_shape_cylinder {
        float radius = 1;
        float zmin   = -1;
        float zmax   = 1;
        float phimax = 360;
    };
    struct pbrt_shape_hyperboloid {
        vec3f p1     = {0, 0, 0};
        vec3f p2     = {1, 1, 1};
        float phimax = 360;
    };
    struct pbrt_shape_paraboloid {
        float radius = 1;
        float zmin   = 0;
        float zmax   = 1;
        float phimax = 360;
    };
    struct pbrt_shape_heightfield {
        int           nu = 0;
        int           nv = 0;
        vector<float> Pz = {};
    };
    if (type == "trianglemesh") {
        auto tvalue = pbrt_shape_trianglemesh{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "indices") {
                parse_param(streams, ptype, tvalue.indices);
            } else if (pname == "P") {
                parse_param(streams, ptype, tvalue.P);
            } else if (pname == "N") {
                parse_param(streams, ptype, tvalue.N);
            } else if (pname == "S") {
                parse_param(streams, ptype, tvalue.S);
            } else if (pname == "uv") {
                parse_param(streams, ptype, tvalue.uv);
            } else if (pname == "alpha") {
                parse_param(streams, ptype, tvalue.alpha);
            } else if (pname == "shadowalpha") {
                parse_param(streams, ptype, tvalue.shadowalpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Shape " + type);
    }
}

// Parse AreaLightSource
void parse_pbrt_arealight(
    vector<pbrt_token_stream>& streams, pbrt_arealight& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
    if (type == "diffuse") {
        auto tvalue = pbrt_arealight_diffuse{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "L") {
                parse_param(streams, ptype, tvalue.L);
            } else if (pname == "twosided") {
                parse_param(streams, ptype, tvalue.twosided);
            } else if (pname == "samples") {
                parse_param(streams, ptype, tvalue.samples);
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
        } else if (cmd == "AttributeBegin") {
        } else if (cmd == "AttributeEnd") {
        } else if (cmd == "Transform") {
            auto xf = identity_mat4f;
            parse_param(streams, xf);
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
        } else if (cmd == "Texture") {
            auto name = ""s;
            parse_value(streams, name);
            auto value = pbrt_texture{};
            parse_pbrt_texture(streams, value);
        } else if (cmd == "Material") {
            static auto material_id = 0;
            auto name  = "unnamed_material_" + std::to_string(material_id++);
            auto value = pbrt_material{};
            parse_pbrt_material(streams, value, false);
        } else if (cmd == "MakeNamedMaterial") {
            auto name = ""s;
            parse_value(streams, name);
            auto value = pbrt_material{};
            parse_pbrt_material(streams, value, true);
        } else if (cmd == "NamedMaterial") {
            auto name = ""s;
            parse_value(streams, name);
        } else if (cmd == "Shape") {
            auto value = pbrt_shape{};
            parse_pbrt_shape(streams, value);
        } else if (cmd == "AreaLightSource") {
            auto value = pbrt_arealight{};
            parse_pbrt_arealight(streams, value);
        } else {
            throw pbrtio_error("unknown command " + cmd);
        }
    }
}

}  // namespace yocto
