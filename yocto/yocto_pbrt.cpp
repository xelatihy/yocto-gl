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
#include "yocto_image.h"

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
    string_view saved;
};

// skip white space or comment
static inline void skip_whitespace_or_comment(pbrt_token_stream& stream) {
    auto& str = stream.str;
    if (str.empty()) return;
    while (!str.empty() && (std::isspace(str.front()) || str.front() == '#' ||
                               str.front() == ',')) {
        if (str.front() == '#') {
            auto pos = str.find('\n');
            if (pos != string_view::npos) {
                str.remove_prefix(pos);
            } else {
                str.remove_prefix(str.length());
            }
        } else {
            auto pos = str.find_first_not_of(" \t\n\r,");
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
template <typename T>
static inline void parse_value(pbrt_token_stream& stream, T& value,
    unordered_map<string, T>& value_names) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("expected enum value");
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
    return parse_value(stream, value, value_names);
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
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_imagemap::wrap_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_imagemap::wrap_type>{
            {"repeat", pbrt_texture_imagemap::wrap_type::repeat},
            {"clamp", pbrt_texture_imagemap::wrap_type::clamp},
            {"black", pbrt_texture_imagemap::wrap_type::black},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_shape_curve::basis_t& value) {
    static auto value_names = unordered_map<string, pbrt_shape_curve::basis_t>{
        {"bezier", pbrt_shape_curve::basis_t::bezier},
        {"bspline", pbrt_shape_curve::basis_t::bspline},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_shape_curve::type_t& value) {
    static auto value_names = unordered_map<string, pbrt_shape_curve::type_t>{
        {"flat", pbrt_shape_curve::type_t::flat},
        {"cylinder", pbrt_shape_curve::type_t::cylinder},
        {"ribbon", pbrt_shape_curve::type_t::ribbon},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_path::lightsamplestrategy_t&  value) {
    static auto value_names =
        unordered_map<string, pbrt_integrator_path::lightsamplestrategy_t>{
            {"power", pbrt_integrator_path::lightsamplestrategy_t::power},
            {"spatial", pbrt_integrator_path::lightsamplestrategy_t::spatial},
            {"uniform", pbrt_integrator_path::lightsamplestrategy_t::uniform},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(pbrt_token_stream&   stream,
    pbrt_integrator_volpath::lightsamplestrategy_t& value) {
    return parse_value(
        stream, (pbrt_integrator_path::lightsamplestrategy_t&)value);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_bdpt::lightsamplestrategy_t&  value) {
    return parse_value(
        stream, (pbrt_integrator_path::lightsamplestrategy_t&)value);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_directlighting::strategy_t&   value) {
    static auto value_names =
        unordered_map<string, pbrt_integrator_directlighting::strategy_t>{
            {"all", pbrt_integrator_directlighting::strategy_t::all},
            {"one", pbrt_integrator_directlighting::strategy_t::one},
        };
    return parse_value(stream, value, value_names);
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
template <typename T, int N>
static inline void parse_value(
    pbrt_token_stream& stream, spectrum<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}

// Check if empty
static inline bool is_empty(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return stream.str.empty();
}

// Check if empty
static inline bool is_string(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '"';
}
static inline bool is_open_bracket(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '[';
}
static inline bool is_close_bracket(pbrt_token_stream& stream) {
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
        return type == "bool";
    } else if constexpr (std::is_same<T, string>::value) {
        return type == "string";
    } else if constexpr (std::is_same<T, vec2f>::value) {
        return type == "point2" || type == "vector2" || type == "float";
    } else if constexpr (std::is_same<T, vec3f>::value) {
        return type == "point3" || type == "vector3" || type == "normal3" ||
               type == "point" || type == "vector" || type == "normal" ||
               type == "float";
    } else if constexpr (std::is_same<T, spectrum3f>::value) {
        return type == "rgb" || type == "spectrum" || type == "blackbody";
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
    pbrt_token_stream& stream, const string& type, spectrum<T, 3>& value) {
    if (type == "rgb") {
        parse_param(stream, value);
    } else if (type == "color") {
        parse_param(stream, value);
    } else if (type == "blackbody") {
        auto blackbody = zero2f;
        parse_param(stream, blackbody);
        (vec3f&)value = blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "spectrum") {
        printf("spectrum  not well supported\n");
        value = {1, 0, 0};
    } else {
        throw pbrtio_error("unsupported spectrum type");
    }
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, vector<T>& value) {
    if (!is_type_compatible<T>(type)) {
        throw pbrtio_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, pbrt_textured<T>& value) {
    if (type == "texture") {
        parse_param(stream, value.texture);
    } else {
        parse_param(stream, value.value);
    }
}

static inline void skip_value(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.front() == '"') {
        str.remove_prefix(1);
        str.remove_prefix(str.find('"') + 1);
    } else {
        str.remove_prefix(str.find_first_of(" \n\t\r],\""));
    }
    skip_whitespace_or_comment(stream);
}

static inline void skip_param(pbrt_token_stream& stream) {
    if (is_open_bracket(stream)) {
        skip_open_bracket(stream);
        while (!is_close_bracket(stream)) skip_value(stream);
        skip_close_bracket(stream);
    } else {
        skip_value(stream);
    }
}

static inline void save_stream_position(pbrt_token_stream& stream) {
    stream.saved = stream.str;
}
static inline void restore_stream_position(pbrt_token_stream& stream) {
    stream.str = stream.saved;
}

// Load a token stream
static inline void load_token_stream(
    const string& filename, pbrt_token_stream& stream) {
    load_text(filename, stream.buffer);
    stream.str = stream.buffer;
}

// operations on token stacks
static inline void init_token_streams(vector<pbrt_token_stream>& streams) {
    streams.reserve(100);
}
static inline void load_token_stream(
    const string& filename, vector<pbrt_token_stream>& streams) {
    streams.emplace_back();
    load_token_stream(filename, streams.back());
}

// Skip whitespace
static inline void skip_whitespace_or_comment(
    vector<pbrt_token_stream>& streams) {
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
static inline bool is_empty(vector<pbrt_token_stream>& streams) {
    skip_whitespace_or_comment(streams);
    return streams.empty();
}
static inline bool is_param(vector<pbrt_token_stream>& streams) {
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

static inline void save_stream_position(vector<pbrt_token_stream>& streams) {
    save_stream_position(streams.back());
}
static inline void restore_stream_position(vector<pbrt_token_stream>& streams) {
    restore_stream_position(streams.back());
}
static inline void skip_param(vector<pbrt_token_stream>& streams) {
    skip_param(streams.back());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Parse Integrator
static void parse_pbrt_integrator(
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
            } else if (pname == "lightsamplestrategy") {
                parse_param(streams, ptype, tvalue.lightsamplestrategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
            // parse_optional_param(streams, "lightsamplestrategy",
            // tvalue.lightsamplestrategy); // TODO: enums
        }
        value = tvalue;
    } else if (type == "volpath") {
        auto tvalue = pbrt_integrator_volpath{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "rrthreshold") {
                parse_param(streams, ptype, tvalue.rrthreshold);
            } else if (pname == "lightsamplestrategy") {
                parse_param(streams, ptype, tvalue.lightsamplestrategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "directlighting") {
        auto tvalue = pbrt_integrator_directlighting{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "strategy") {
                parse_param(streams, ptype, tvalue.strategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bdpt") {
        auto tvalue = pbrt_integrator_bdpt{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "lightsamplestrategy") {
                parse_param(streams, ptype, tvalue.lightsamplestrategy);
            } else if (pname == "visualizestrategies") {
                parse_param(streams, ptype, tvalue.visualizestrategies);
            } else if (pname == "visualizeweights") {
                parse_param(streams, ptype, tvalue.visualizeweights);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mlt") {
        auto tvalue = pbrt_integrator_mlt{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "bootstrapsamples") {
                parse_param(streams, ptype, tvalue.bootstrapsamples);
            } else if (pname == "chains") {
                parse_param(streams, ptype, tvalue.chains);
            } else if (pname == "mutationsperpixel") {
                parse_param(streams, ptype, tvalue.mutationsperpixel);
            } else if (pname == "largestepprobability") {
                parse_param(streams, ptype, tvalue.largestepprobability);
            } else if (pname == "sigma") {
                parse_param(streams, ptype, tvalue.sigma);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sppm") {
        auto tvalue = pbrt_integrator_sppm{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else if (pname == "iterations") {
                parse_param(streams, ptype, tvalue.iterations);
            } else if (pname == "photonsperiteration") {
                parse_param(streams, ptype, tvalue.photonsperiteration);
            } else if (pname == "imagewritefrequency") {
                parse_param(streams, ptype, tvalue.imagewritefrequency);
            } else if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "whitted") {
        auto tvalue = pbrt_integrator_whitted{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(streams, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(streams, ptype, tvalue.pixelbounds);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Integrator " + type);
    }
}

// Parse Sampler
static void parse_pbrt_sampler(
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
static void parse_pbrt_filter(
    vector<pbrt_token_stream>& streams, pbrt_filter& value) {
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
static void parse_pbrt_film(
    vector<pbrt_token_stream>& streams, pbrt_film& value) {
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
static void parse_pbrt_camera(
    vector<pbrt_token_stream>& streams, pbrt_camera& value) {
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
static void parse_pbrt_texture(
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
        throw pbrtio_error("unknown Texture " + type);
    }
}

// Get typename
static void get_pbrt_typeparam(
    vector<pbrt_token_stream>& streams, string& value) {
    save_stream_position(streams);
    value      = "";
    auto pname = ""s, ptype = ""s;
    while (is_param(streams) && value == "") {
        parse_nametype(streams, pname, ptype);
        if (pname == "type") {
            parse_param(streams, ptype, value);
        } else {
            skip_param(streams);
        }
    }
    if (value == "") throw pbrtio_error("type not found");
    restore_stream_position(streams);
}

// Parse Material
static void parse_pbrt_material(
    vector<pbrt_token_stream>& streams, pbrt_material& value, bool named) {
    auto type = ""s;
    if (!named) {
        parse_value(streams, type);
    } else {
        get_pbrt_typeparam(streams, type);
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
                parse_param(streams, ptype, tvalue.eta);
            } else if (pname == "beta_m") {
                parse_param(streams, ptype, tvalue.beta_m);
            } else if (pname == "beta_n") {
                parse_param(streams, ptype, tvalue.beta_n);
            } else if (pname == "alpha") {
                parse_param(streams, ptype, tvalue.alpha);
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
            } else if (pname == "index") {
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
            } else if (pname == "bumpmap") {
                parse_param(streams, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(streams, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
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
static void parse_pbrt_shape(
    vector<pbrt_token_stream>& streams, pbrt_shape& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
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
    } else if (type == "plymesh") {
        auto tvalue = pbrt_shape_plymesh{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "filename") {
                parse_param(streams, ptype, tvalue.filename);
            } else if (pname == "alpha") {
                parse_param(streams, ptype, tvalue.alpha);
            } else if (pname == "shadowalpha") {
                parse_param(streams, ptype, tvalue.shadowalpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "curve") {
        auto tvalue = pbrt_shape_curve{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "P") {
                parse_param(streams, ptype, tvalue.P);
            } else if (pname == "N") {
                parse_param(streams, ptype, tvalue.N);
            } else if (pname == "basis") {
                parse_param(streams, ptype, tvalue.basis);
            } else if (pname == "degree") {
                parse_param(streams, ptype, tvalue.degree);
            } else if (pname == "type") {
                parse_param(streams, ptype, tvalue.type);
            } else if (pname == "width") {
                auto width = 1.0f;
                parse_param(streams, ptype, width);
                tvalue.width0 = width;
                tvalue.width1 = width;
            } else if (pname == "width0") {
                parse_param(streams, ptype, tvalue.width0);
            } else if (pname == "width1") {
                parse_param(streams, ptype, tvalue.width1);
            } else if (pname == "splitdepth") {
                parse_param(streams, ptype, tvalue.splitdepth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "loopsubdiv") {
        auto tvalue = pbrt_shape_loopsubdiv{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "indices") {
                parse_param(streams, ptype, tvalue.indices);
            } else if (pname == "P") {
                parse_param(streams, ptype, tvalue.P);
            } else if (pname == "levels") {
                parse_param(streams, ptype, tvalue.levels);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "nurbs") {
        auto tvalue = pbrt_shape_nurbs{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "nu") {
                parse_param(streams, ptype, tvalue.nu);
            } else if (pname == "nv") {
                parse_param(streams, ptype, tvalue.nv);
            } else if (pname == "uknots") {
                parse_param(streams, ptype, tvalue.uknots);
            } else if (pname == "vknots") {
                parse_param(streams, ptype, tvalue.vknots);
            } else if (pname == "u0") {
                parse_param(streams, ptype, tvalue.u0);
            } else if (pname == "v0") {
                parse_param(streams, ptype, tvalue.v0);
            } else if (pname == "u1") {
                parse_param(streams, ptype, tvalue.u1);
            } else if (pname == "v1") {
                parse_param(streams, ptype, tvalue.v1);
            } else if (pname == "P") {
                parse_param(streams, ptype, tvalue.P);
            } else if (pname == "Pw") {
                parse_param(streams, ptype, tvalue.Pw);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sphere") {
        auto tvalue = pbrt_shape_sphere{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(streams, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(streams, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disk") {
        auto tvalue = pbrt_shape_disk{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else if (pname == "height") {
                parse_param(streams, ptype, tvalue.height);
            } else if (pname == "innerradius") {
                parse_param(streams, ptype, tvalue.innerradius);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cone") {
        auto tvalue = pbrt_shape_cone{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else if (pname == "height") {
                parse_param(streams, ptype, tvalue.height);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cylinder") {
        auto tvalue = pbrt_shape_cylinder{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(streams, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(streams, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hyperboloid") {
        auto tvalue = pbrt_shape_hyperboloid{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "p1") {
                parse_param(streams, ptype, tvalue.p1);
            } else if (pname == "p2") {
                parse_param(streams, ptype, tvalue.p2);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "paraboloid") {
        auto tvalue = pbrt_shape_paraboloid{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "radius") {
                parse_param(streams, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(streams, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(streams, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(streams, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heightfield") {
        auto tvalue = pbrt_shape_heightfield{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "nu") {
                parse_param(streams, ptype, tvalue.nu);
            } else if (pname == "nv") {
                parse_param(streams, ptype, tvalue.nv);
            } else if (pname == "Pz") {
                parse_param(streams, ptype, tvalue.Pz);
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
static void parse_pbrt_arealight(
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

// Parse LightSource
static void parse_pbrt_light(
    vector<pbrt_token_stream>& streams, pbrt_light& value) {
    auto type = ""s;
    parse_value(streams, type);
    auto pname = ""s, ptype = ""s;
    if (type == "distant") {
        auto tvalue = pbrt_light_distant{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(streams, ptype, tvalue.L);
            } else if (pname == "from") {
                parse_param(streams, ptype, tvalue.from);
            } else if (pname == "to") {
                parse_param(streams, ptype, tvalue.to);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "goniometric") {
        auto tvalue = pbrt_light_goniometric{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(streams, ptype, tvalue.I);
            } else if (pname == "mapname") {
                parse_param(streams, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "infinite") {
        auto tvalue = pbrt_light_infinite{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(streams, ptype, tvalue.L);
            } else if (pname == "samples") {
                parse_param(streams, ptype, tvalue.samples);
            } else if (pname == "mapname") {
                parse_param(streams, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "distant") {
        auto tvalue = pbrt_light_distant{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(streams, ptype, tvalue.L);
            } else if (pname == "from") {
                parse_param(streams, ptype, tvalue.from);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "projection") {
        auto tvalue = pbrt_light_projection{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(streams, ptype, tvalue.I);
            } else if (pname == "fov") {
                parse_param(streams, ptype, tvalue.fov);
            } else if (pname == "mapname") {
                parse_param(streams, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "spot") {
        auto tvalue = pbrt_light_spot{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(streams, ptype, tvalue.I);
            } else if (pname == "from") {
                parse_param(streams, ptype, tvalue.from);
            } else if (pname == "to") {
                parse_param(streams, ptype, tvalue.to);
            } else if (pname == "coneangle") {
                parse_param(streams, ptype, tvalue.coneangle);
            } else if (pname == "conedeltaangle") {
                parse_param(streams, ptype, tvalue.conedeltaangle);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Medium
static void parse_pbrt_medium(
    vector<pbrt_token_stream>& streams, pbrt_medium& value, bool named) {
    auto type = ""s;
    if (!named) {
        parse_value(streams, type);
    } else {
        get_pbrt_typeparam(streams, type);
    }
    auto pname = ""s, ptype = ""s;
    if (type == "distant") {
        auto tvalue = pbrt_medium_homogeneous{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "sigma_a") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "sigma_s") {
                parse_param(streams, ptype, tvalue.sigma_s);
            } else if (pname == "preset") {
                parse_param(streams, ptype, tvalue.preset);
            } else if (pname == "g") {
                parse_param(streams, ptype, tvalue.g);
            } else if (pname == "scale") {
                parse_param(streams, ptype, tvalue.scale);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heterogeneous") {
        auto tvalue = pbrt_medium_heterogeneous{};
        while (is_param(streams)) {
            parse_nametype(streams, pname, ptype);
            if (pname == "sigma_a") {
                parse_param(streams, ptype, tvalue.scale);
            } else if (pname == "sigma_s") {
                parse_param(streams, ptype, tvalue.sigma_s);
            } else if (pname == "preset") {
                parse_param(streams, ptype, tvalue.preset);
            } else if (pname == "g") {
                parse_param(streams, ptype, tvalue.g);
            } else if (pname == "p0") {
                parse_param(streams, ptype, tvalue.p0);
            } else if (pname == "p1") {
                parse_param(streams, ptype, tvalue.p1);
            } else if (pname == "nx") {
                parse_param(streams, ptype, tvalue.nx);
            } else if (pname == "ny") {
                parse_param(streams, ptype, tvalue.ny);
            } else if (pname == "nz") {
                parse_param(streams, ptype, tvalue.nz);
            } else if (pname == "density") {
                parse_param(streams, ptype, tvalue.density);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        struct pbrt_heterogeneous_medium {
            vector<float> density = {};
        };
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

    // parsing stack
    auto stack  = vector<pbrt_context>{{}};
    auto object = pbrt_object{};

    // parse command by command
    auto cmd = ""s;
    while (!is_empty(streams)) {
        // get command
        parse_command(streams, cmd);
        if (cmd == "WorldBegin") {
            stack.push_back({});
        } else if (cmd == "WorldEnd") {
            stack.pop_back();
            if (stack.size() != 1) throw pbrtio_error("bad stack");
        } else if (cmd == "AttributeBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "AttributeEnd") {
            stack.pop_back();
        } else if (cmd == "TransformBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "TransformEnd") {
            stack.pop_back();
        } else if (cmd == "ObjectBegin") {
            parse_value(streams, object.name);
            if (cb.begin_object) cb.begin_object(object, stack.back());
        } else if (cmd == "ObjectEnd") {
            if (cb.end_object) cb.end_object(object, stack.back());
            object = {};
        } else if (cmd == "ObjectInstance") {
            auto value = pbrt_object{};
            parse_value(streams, value.name);
            if (cb.object_instance) cb.object_instance(value, stack.back());
        } else if (cmd == "Transform") {
            auto xf = identity_mat4f;
            parse_param(streams, xf);
            stack.back().frame = (affine3f)xf;
        } else if (cmd == "ConcatTransform") {
            auto xf = identity_mat4f;
            parse_param(streams, xf);
            stack.back().frame = stack.back().frame * (affine3f)xf;
        } else if (cmd == "Scale") {
            auto v = zero3f;
            parse_param(streams, v);
            stack.back().frame = stack.back().frame *
                                 (affine3f)make_scaling_frame(v);
        } else if (cmd == "Translate") {
            auto v = zero3f;
            parse_param(streams, v);
            stack.back().frame = stack.back().frame *
                                 (affine3f)make_translation_frame(v);
        } else if (cmd == "Rotate") {
            auto v = zero4f;
            parse_param(streams, v);
            stack.back().frame = stack.back().frame *
                                 (affine3f)(affine3f)make_rotation_frame(
                                     vec3f{v.y, v.z, v.w}, radians(v.x));
        } else if (cmd == "LookAt") {
            auto from = zero3f, to = zero3f, up = zero3f;
            parse_param(streams, from);
            parse_param(streams, to);
            parse_param(streams, up);
            // from pbrt parser
            auto frame = make_lookat_frame(from, to, up, true);
            // frame.z = normalize(to-from);
            // frame.x = normalize(cross(frame.z,up));
            // frame.y = cross(frame.x,frame.z);
            // frame.o    = from;
            stack.back().frame = stack.back().frame * (affine3f)inverse(frame);
            // stack.back().focus = length(m.x - m.y);
        } else if (cmd == "ReverseOrientation") {
            stack.back().reverse = !stack.back().reverse;
        } else if (cmd == "Integrator") {
            auto value = pbrt_integrator{};
            parse_pbrt_integrator(streams, value);
            if (cb.integrator) cb.integrator(value, stack.back());
        } else if (cmd == "Sampler") {
            auto value = pbrt_sampler{};
            parse_pbrt_sampler(streams, value);
            if (cb.sampler) cb.sampler(value, stack.back());
        } else if (cmd == "PixelFilter") {
            auto value = pbrt_filter{};
            parse_pbrt_filter(streams, value);
            if (cb.filter) cb.filter(value, stack.back());
        } else if (cmd == "Film") {
            auto value = pbrt_film{};
            parse_pbrt_film(streams, value);
            if (cb.film) cb.film(value, stack.back());
        } else if (cmd == "Camera") {
            auto value = pbrt_camera{};
            parse_pbrt_camera(streams, value);
            if (cb.camera) cb.camera(value, stack.back());
        } else if (cmd == "Texture") {
            auto name = ""s, type = ""s;
            parse_value(streams, name);
            parse_value(streams, type);
            auto value = pbrt_texture{};
            parse_pbrt_texture(streams, value);
            if (cb.texture) cb.texture(value, name, stack.back());
        } else if (cmd == "Material") {
            static auto material_id = 0;
            auto name  = "unnamed_material_" + std::to_string(material_id++);
            auto value = pbrt_material{};
            parse_pbrt_material(streams, value, false);
            stack.back().material = name;
            if (cb.material) cb.material(value, name, stack.back());
        } else if (cmd == "MakeNamedMaterial") {
            auto name = ""s;
            parse_value(streams, name);
            auto value = pbrt_material{};
            parse_pbrt_material(streams, value, true);
            if (cb.material) cb.material(value, name, stack.back());
        } else if (cmd == "NamedMaterial") {
            auto name = ""s;
            parse_value(streams, name);
            stack.back().material = name;
        } else if (cmd == "Shape") {
            auto value = pbrt_shape{};
            parse_pbrt_shape(streams, value);
            if (cb.shape) cb.shape(value, stack.back());
        } else if (cmd == "AreaLightSource") {
            static auto material_id = 0;
            auto name  = "unnamed_arealight_" + std::to_string(material_id++);
            auto value = pbrt_arealight{};
            parse_pbrt_arealight(streams, value);
            stack.back().arealight = name;
            if (cb.arealight) cb.arealight(value, name, stack.back());
        } else if (cmd == "LightSource") {
            auto value = pbrt_light{};
            parse_pbrt_light(streams, value);
            if (cb.light) cb.light(value, stack.back());
        } else if (cmd == "MakeNamedMedium") {
            auto name = ""s;
            parse_value(streams, name);
            auto value = pbrt_medium{};
            parse_pbrt_medium(streams, value, true);
            if (cb.medium) cb.medium(value, name, stack.back());
        } else if (cmd == "MediumInterface") {
            auto interior = ""s, exterior = ""s;
            parse_value(streams, interior);
            parse_value(streams, exterior);
            stack.back().medium_interior = interior;
            stack.back().medium_exterior = exterior;
        } else {
            throw pbrtio_error("unknown command " + cmd);
        }
    }
}

}  // namespace yocto
