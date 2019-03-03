//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
#include "ysceneui.h"

#include "ext/CLI11.hpp"

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};

void exit_error(const string& msg) {
    printf("%s\n", msg.c_str());
    exit(1);
}

struct drawgl_shape {
    opengl_array_buffer  positions_buffer     = {};
    opengl_array_buffer  normals_buffer       = {};
    opengl_array_buffer  texturecoords_buffer = {};
    opengl_array_buffer  colors_buffer        = {};
    opengl_array_buffer  tangentspaces_buffer = {};
    opengl_elementbuffer points_buffer        = {};
    opengl_elementbuffer lines_buffer         = {};
    opengl_elementbuffer triangles_buffer     = {};
    opengl_elementbuffer quads_buffer         = {};
};

struct drawgl_state {
    opengl_program         program  = {};
    vector<drawgl_shape>   shapes   = {};
    vector<opengl_texture> textures = {};
};

struct drawgl_lights {
    vector<vec3f> positions = {};
    vector<vec3f> emission  = {};
    vector<int>   types     = {};

    bool empty() const { return positions.empty(); }
};

drawgl_lights make_drawgl_lights(const yocto_scene& scene) {
    auto lights = drawgl_lights{};
    for (auto& instance : scene.instances) {
        if (instance.shape < 0) continue;
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[shape.material];
        if (material.emission == zero3f) continue;
        if (lights.positions.size() >= 16) break;
        auto bbox = compute_shape_bounds(shape);
        auto pos  = (bbox.max + bbox.min) / 2;
        auto area = 0.0f;
        if (!shape.triangles.empty()) {
            for (auto t : shape.triangles)
                area += triangle_area(shape.positions[t.x],
                    shape.positions[t.y], shape.positions[t.z]);
        } else if (!shape.quads.empty()) {
            for (auto q : shape.quads)
                area += quad_area(shape.positions[q.x], shape.positions[q.y],
                    shape.positions[q.z], shape.positions[q.w]);
        } else if (!shape.lines.empty()) {
            for (auto l : shape.lines)
                area += line_length(shape.positions[l.x], shape.positions[l.y]);
        } else {
            area += shape.positions.size();
        }
        auto ke = material.emission * area;
        lights.positions.push_back(transform_point(instance.frame, pos));
        lights.emission.push_back(ke);
        lights.types.push_back(0);
    }
    return lights;
}

// Draw options
struct drawgl_options {
    int   camera_id    = 0;
    int   image_width  = 1280;
    int   image_height = 720;
    bool  wireframe    = false;
    bool  edges        = false;
    float edge_offset  = 0.01f;
    bool  eyelight     = false;
    float exposure     = 0;
    float gamma        = 2.2f;
    vec3f ambient      = {0, 0, 0};
    bool  double_sided = false;
    float near_plane   = 0.01f;
    float far_plane    = 10000.0f;
};

// Application state
struct app_state {
    // loading parameters
    string filename     = "scene.json";
    string imfilename   = "out.png";
    string outfilename  = "scene.json";
    bool   double_sided = false;

    // options
    load_scene_options load_options = {};
    build_bvh_options  bvh_options  = {};
    drawgl_options     draw_options = {};

    // scene
    yocto_scene scene = {};

    // rendering state
    drawgl_state  state  = {};
    drawgl_lights lights = {};

    // view image
    bool                          widgets_open   = false;
    bool                          navigation_fps = false;
    pair<type_index, int>         selection      = {typeid(void), -1};
    vector<pair<type_index, int>> update_list;
    float                         time       = 0;
    string                        anim_group = "";
    vec2f                         time_range = zero2f;
    bool                          animate    = false;

    // app status
    bool   load_done = false, load_running = false;
    string status = "";
};

bool load_scene_sync(app_state& app) {
    // scene loading
    app.status = "loading scene";
    try {
        load_scene(app.filename, app.scene, app.load_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
        return false;
    }

    // tesselate
    app.status = "tesselating surfaces";
    tesselate_shapes(app.scene);

    // add components
    add_missing_cameras(app.scene);
    add_missing_names(app.scene);
    print_validation_errors(app.scene);

    // init renderer
    app.status = "initializing lights";
    app.lights = make_drawgl_lights(app.scene);

    // fix renderer type if no lights
    if (app.lights.empty() && !app.draw_options.eyelight) {
        printf("no lights presents, switching to eyelight shader\n");
        app.draw_options.eyelight = true;
    }

    // animation
    app.time_range = compute_animation_range(app.scene);
    app.time       = app.time_range.x;

    // set flags
    app.load_done    = true;
    app.load_running = false;
    app.status       = "loading done";

    // done
    return false;
}

void load_scene_async(app_state& app) {
    if (app.load_running) {
        throw runtime_error("already loading");
    }
    app.load_done    = false;
    app.load_running = true;
    app.status       = "uninitialized";
    app.scene        = {};
    app.lights       = {};
    app.state        = {};
    auto load_thread = thread([&app]() { load_scene_sync(app); });
    load_thread.detach();
}

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

static const char* vertex =
    R"(
        #version 330

        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;           // shape transform
        uniform float shape_normal_offset;           // shape normal offset

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        out vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
        out vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
        out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
        out vec4 color;                 // [to fragment shader] vertex color
        out vec4 tangsp;                // [to fragment shader] vertex tangent space

        // main function
        void main() {
            // copy values
            pos = vert_pos;
            norm = vert_norm;
            tangsp = vert_tangsp;

            // normal offset
            if(shape_normal_offset != 0) {
                pos += shape_normal_offset * norm;
            }

            // world projection
            pos = (shape_xform * vec4(pos,1)).xyz;
            norm = (shape_xform * vec4(norm,0)).xyz;
            tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
        }
        )";

static const char* fragment =
    R"(
        #version 330

        float pif = 3.14159265;

        uniform bool eyelight;         // eyelight shading
        uniform vec3 lamb;             // ambient light
        uniform int lnum;              // number of lights
        uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
        uniform vec3 lpos[16];         // light positions
        uniform vec3 lke[16];          // light intensities

        void evaluate_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(ltype[lid] == 0) {
                // compute point light color at pos
                cl = lke[lid] / pow(length(lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lpos[lid]-pos);
            }
            else if(ltype[lid] == 1) {
                // compute light color
                cl = lke[lid];
                // compute light direction
                wi = normalize(lpos[lid]);
            }
        }

        vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
            vec3 n, vec3 wi, vec3 wo) {
            if(etype == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(rs*rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(etype == 1) {
                return ((1+dot(wo,wi))/2) * kd/pif;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pif;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pif)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pif;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pif)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        uniform int elem_type;
        uniform bool elem_faceted;
        uniform vec4 highlight;   // highlighted color

        uniform int mat_type;          // material type
        uniform vec3 mat_ke;           // material ke
        uniform vec3 mat_kd;           // material kd
        uniform vec3 mat_ks;           // material ks
        uniform float mat_rs;          // material rs
        uniform float mat_op;          // material op

        uniform bool mat_ke_txt_on;    // material ke texture on
        uniform sampler2D mat_ke_txt;  // material ke texture
        uniform bool mat_kd_txt_on;    // material kd texture on
        uniform sampler2D mat_kd_txt;  // material kd texture
        uniform bool mat_ks_txt_on;    // material ks texture on
        uniform sampler2D mat_ks_txt;  // material ks texture
        uniform bool mat_rs_txt_on;    // material rs texture on
        uniform sampler2D mat_rs_txt;  // material rs texture
        uniform bool mat_op_txt_on;    // material op texture on
        uniform sampler2D mat_op_txt;  // material op texture

        uniform bool mat_norm_txt_on;    // material norm texture on
        uniform sampler2D mat_norm_txt;  // material norm texture

        uniform bool mat_double_sided;   // double sided rendering

        uniform mat4 shape_xform;           // shape transform

        bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                           out vec3 kd, out vec3 ks, out float rs, out float op) {
            if(mat_type == 0) {
                ke = mat_ke;
                kd = vec3(0,0,0);
                ks = vec3(0,0,0);
                op = 1;
                return false;
            }

            ke = color.xyz * mat_ke;
            kd = color.xyz * mat_kd;
            ks = color.xyz * mat_ks;
            rs = mat_rs;
            op = color.w * mat_op;

            vec4 ke_txt = (mat_ke_txt_on) ? texture(mat_ke_txt,texcoord) : vec4(1,1,1,1);
            vec4 kd_txt = (mat_kd_txt_on) ? texture(mat_kd_txt,texcoord) : vec4(1,1,1,1);
            vec4 ks_txt = (mat_ks_txt_on) ? texture(mat_ks_txt,texcoord) : vec4(1,1,1,1);
            vec4 rs_txt = (mat_rs_txt_on) ? texture(mat_rs_txt,texcoord) : vec4(1,1,1,1);
            vec4 op_txt = (mat_op_txt_on) ? texture(mat_op_txt,texcoord) : vec4(1,1,1,1);

            // get material color from textures and adjust values
            if(mat_type == 1) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt.y;
                rs = rs*rs;
                op *= op_txt.x * kd_txt.w;
            } else if(mat_type == 2) {
                ke *= ke_txt.xyz;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(mat_type == 3) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                float gs = (1 - rs) * ks_txt.w;
                rs = 1 - gs;
                rs = rs*rs;
                op *= kd_txt.w;
            }

            return true;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!mat_norm_txt_on) return norm;
            vec3 tangu = normalize((shape_xform * vec4(normalize(tangsp.xyz),0)).xyz);
            vec3 tangv = normalize(cross(norm, tangu));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 texture = 2 * texture(mat_norm_txt,texcoord).xyz - 1;
            texture.y = -texture.y;
            return normalize( tangu * texture.x + tangv * texture.y + norm * texture.z );
        }

        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        uniform vec3 cam_pos;          // camera position
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        uniform float exposure; 
        uniform float gamma;

        out vec4 frag_color;      

        vec3 triangle_normal(vec3 pos) {
            vec3 fdx = dFdx(pos); 
            vec3 fdy = dFdy(pos); 
            return normalize((shape_xform * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
        }

        // main
        void main() {
            // view vector
            vec3 wo = normalize(cam_pos - pos);

            // prepare normals
            vec3 n;
            if(elem_faceted) {
                n = triangle_normal(pos);
            } else {
                n = normalize(norm);
            }

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(mat_double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
            bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

            // exit if needed
            if(brdf_op < 0.005) discard;

            // check const color
            if(elem_type == 0) {
                frag_color = vec4(brdf_ke,brdf_op);
                return;
            }

            // emission
            vec3 c = brdf_ke;

            // check early exit
            if(brdf_kd != vec3(0,0,0) || brdf_ks != vec3(0,0,0)) {
                // eyelight shading
                if(eyelight) {
                    vec3 wi = wo;
                    c += pif * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lamb * brdf_kd;
                    // foreach light
                    for(int lid = 0; lid < lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        evaluate_light(lid, pos, cl, wi);
                        c += cl * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                    }
                }
            }

            // final color correction
            c = pow(c * pow(2,exposure), vec3(1/gamma));

            // highlighting
            if(highlight.w > 0) {
                if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
                    c = highlight.xyz * highlight.w + c * (1-highlight.w);
            }

            // output final color by setting gl_FragColor
            frag_color = vec4(c,brdf_op);
        }
        )";
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

// Draw a shape
void draw_glinstance(drawgl_state& state, const yocto_scene& scene,
    const yocto_instance& instance, bool highlighted,
    const drawgl_options& options) {
    auto& shape    = scene.shapes[instance.shape];
    auto& vbos     = state.shapes.at(instance.shape);
    auto& material = scene.materials[shape.material];

    auto xform = frame_to_mat(instance.frame);

    set_opengl_uniform(state.program, "shape_xform", xform);
    set_opengl_uniform(state.program, "shape_normal_offset", 0.0f);
    set_opengl_uniform(
        state.program, "highlight", (highlighted) ? vec4f{1, 1, 0, 1} : zero4f);

    auto mtype = 1;
    if (material.base_metallic) mtype = 2;
    if (material.gltf_textures) mtype = (material.base_metallic) ? 2 : 3;
    set_opengl_uniform(state.program, "mat_type", mtype);
    set_opengl_uniform(state.program, "mat_ke", material.emission);
    set_opengl_uniform(state.program, "mat_kd", material.diffuse);
    set_opengl_uniform(state.program, "mat_ks", material.specular);
    set_opengl_uniform(state.program, "mat_rs", material.roughness);
    set_opengl_uniform(state.program, "mat_op", material.opacity);
    set_opengl_uniform(
        state.program, "mat_double_sided", (int)options.double_sided);
    set_opengl_uniform_texture(state.program, "mat_ke_txt", "mat_ke_txt_on",
        material.emission_texture >= 0
            ? state.textures.at(material.emission_texture)
            : opengl_texture{},
        0);
    set_opengl_uniform_texture(state.program, "mat_kd_txt", "mat_kd_txt_on",
        material.diffuse_texture >= 0
            ? state.textures.at(material.diffuse_texture)
            : opengl_texture{},
        1);
    set_opengl_uniform_texture(state.program, "mat_ks_txt", "mat_ks_txt_on",
        material.specular_texture >= 0
            ? state.textures.at(material.specular_texture)
            : opengl_texture{},
        2);
    set_opengl_uniform_texture(state.program, "mat_rs_txt", "mat_rs_txt_on",
        material.roughness_texture >= 0
            ? state.textures.at(material.roughness_texture)
            : opengl_texture{},
        3);
    set_opengl_uniform_texture(state.program, "mat_norm_txt", "mat_norm_txt_on",
        material.normal_texture >= 0
            ? state.textures.at(material.normal_texture)
            : opengl_texture{},
        5);

    set_opengl_uniform(
        state.program, "elem_faceted", (int)shape.normals.empty());
    set_opengl_vertexattrib(
        state.program, "vert_pos", vbos.positions_buffer, zero3f);
    set_opengl_vertexattrib(
        state.program, "vert_norm", vbos.normals_buffer, zero3f);
    set_opengl_vertexattrib(
        state.program, "vert_texcoord", vbos.texturecoords_buffer, zero2f);
    set_opengl_vertexattrib(
        state.program, "vert_color", vbos.colors_buffer, vec4f{1, 1, 1, 1});
    set_opengl_vertexattrib(state.program, "vert_tangsp",
        vbos.tangentspaces_buffer, vec4f{0, 0, 1, 1});

    if (vbos.points_buffer) {
        set_opengl_uniform(state.program, "elem_type", 1);
        draw_opengl_points(vbos.points_buffer, vbos.points_buffer.num);
    }
    if (vbos.lines_buffer) {
        set_opengl_uniform(state.program, "elem_type", 2);
        draw_opengl_lines(vbos.lines_buffer, vbos.lines_buffer.num);
    }
    if (vbos.triangles_buffer) {
        set_opengl_uniform(state.program, "elem_type", 3);
        draw_opengl_triangles(vbos.triangles_buffer, vbos.triangles_buffer.num);
    }
    if (vbos.quads_buffer) {
        set_opengl_uniform(state.program, "elem_type", 3);
        draw_opengl_triangles(vbos.quads_buffer, vbos.quads_buffer.num);
    }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_opengl_error();
        set_opengl_uniform(state.program, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state.program, "ke"), 0, 0, 0);
        set_opengl_uniform(state.program, "op"), material.op);
        set_opengl_uniform(state.program, "shp_normal_offset"), 0.01f);
        check_opengl_error();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_opengl_error();
    }
#endif
    if (options.edges) throw runtime_error("edges are momentarily disabled");

    // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
}

// Display a scene
void draw_glscene(drawgl_state& state, const yocto_scene& scene,
    const vec2i& viewport_size, const pair<type_index, int>& highlighted,
    const drawgl_options& options) {
    auto& camera      = scene.cameras.at(options.camera_id);
    auto  camera_view = frame_to_mat(inverse(camera.frame));
    auto  camera_proj = make_perspective_mat(get_camera_fovy(camera),
        (float)viewport_size.x / (float)viewport_size.y, options.near_plane,
        options.far_plane);

    bind_opengl_program(state.program);
    set_opengl_uniform(state.program, "cam_pos", camera.frame.o);
    set_opengl_uniform(state.program, "cam_xform_inv", camera_view);
    set_opengl_uniform(state.program, "cam_proj", camera_proj);
    set_opengl_uniform(state.program, "eyelight", (int)options.eyelight);
    set_opengl_uniform(state.program, "exposure", options.exposure);
    set_opengl_uniform(state.program, "gamma", options.gamma);

    if (!options.eyelight) {
        auto lights_pos  = vector<vec3f>();
        auto lights_ke   = vector<vec3f>();
        auto lights_type = vector<int>();
        for (auto& instance : scene.instances) {
            if (instance.shape < 0) continue;
            auto& shape    = scene.shapes[instance.shape];
            auto& material = scene.materials[shape.material];
            if (material.emission == zero3f) continue;
            if (lights_pos.size() >= 16) break;
            auto bbox = compute_shape_bounds(shape);
            auto pos  = (bbox.max + bbox.min) / 2;
            auto area = 0.0f;
            if (!shape.triangles.empty()) {
                for (auto t : shape.triangles)
                    area += triangle_area(shape.positions[t.x],
                        shape.positions[t.y], shape.positions[t.z]);
            } else if (!shape.quads.empty()) {
                for (auto q : shape.quads)
                    area += quad_area(shape.positions[q.x],
                        shape.positions[q.y], shape.positions[q.z],
                        shape.positions[q.w]);
            } else if (!shape.lines.empty()) {
                for (auto l : shape.lines)
                    area += line_length(
                        shape.positions[l.x], shape.positions[l.y]);
            } else {
                area += shape.positions.size();
            }
            auto ke = material.emission * area;
            lights_pos.push_back(transform_point(instance.frame, pos));
            lights_ke.push_back(ke);
            lights_type.push_back(0);
        }
        set_opengl_uniform(state.program, "lamb", zero3f);
        set_opengl_uniform(state.program, "lnum", (int)lights_pos.size());
        for (auto i = 0; i < lights_pos.size(); i++) {
            auto is = std::to_string(i);
            set_opengl_uniform(
                state.program, ("lpos[" + is + "]").c_str(), lights_pos[i]);
            set_opengl_uniform(
                state.program, ("lke[" + is + "]").c_str(), lights_ke[i]);
            set_opengl_uniform(state.program, ("ltype[" + is + "]").c_str(),
                (int)lights_type[i]);
        }
    }

    if (options.wireframe) set_opengl_wireframe(true);
    for (auto instance_id = 0; instance_id < scene.instances.size();
         instance_id++) {
        auto& instance = scene.instances[instance_id];
        // auto& shape     = scene.shapes[instance.shape];
        // auto& material  = scene.materials[shape.material];
        auto highlight = highlighted ==
                         pair<type_index, int>{
                             typeid(yocto_instance), instance_id};
        draw_glinstance(state, scene, instance, highlight, options);
    }

    unbind_opengl_program();
    if (options.wireframe) set_opengl_wireframe(false);
}

void init_drawgl_state(drawgl_state& state, const yocto_scene& scene) {
    // load textures and vbos
    init_opengl_program(state.program, vertex, fragment);
    state.textures.resize(scene.textures.size());
    for (auto texture_id = 0; texture_id < scene.textures.size();
         texture_id++) {
        auto& texture = scene.textures[texture_id];
        if (!texture.hdr_image.empty()) {
            init_opengl_texture(state.textures[texture_id], texture.hdr_image,
                true, true, true);
        } else if (!texture.ldr_image.empty()) {
            init_opengl_texture(state.textures[texture_id], texture.ldr_image,
                !texture.ldr_as_linear, true, true);
        } else {
            throw runtime_error("bad texture");
        }
    }
    state.shapes.resize(scene.shapes.size());
    for (auto shape_id = 0; shape_id < scene.shapes.size(); shape_id++) {
        auto& shape = scene.shapes[shape_id];
        auto  vbos  = drawgl_shape();
        if (shape.quads_positions.empty()) {
            if (!shape.positions.empty())
                init_opengl_array_buffer(
                    vbos.positions_buffer, shape.positions, false);
            if (!shape.normals.empty())
                init_opengl_array_buffer(
                    vbos.normals_buffer, shape.normals, false);
            if (!shape.texturecoords.empty())
                init_opengl_array_buffer(
                    vbos.texturecoords_buffer, shape.texturecoords, false);
            if (!shape.colors.empty())
                init_opengl_array_buffer(
                    vbos.colors_buffer, shape.colors, false);
            if (!shape.tangentspaces.empty())
                init_opengl_array_buffer(
                    vbos.tangentspaces_buffer, shape.tangentspaces, false);
            if (!shape.points.empty())
                init_opengl_elementbuffer(
                    vbos.points_buffer, shape.points, false);
            if (!shape.lines.empty())
                init_opengl_elementbuffer(
                    vbos.lines_buffer, shape.lines, false);
            if (!shape.triangles.empty())
                init_opengl_elementbuffer(
                    vbos.triangles_buffer, shape.triangles, false);
            if (!shape.quads.empty()) {
                auto triangles = vector<vec3i>{};
                convert_quads_to_triangles(triangles, shape.quads);
                init_opengl_elementbuffer(vbos.quads_buffer, triangles, false);
            }
        } else {
            auto quads         = vector<vec4i>{};
            auto positions     = vector<vec3f>{};
            auto normals       = vector<vec3f>{};
            auto texturecoords = vector<vec2f>{};
            convert_facevarying(quads, positions, normals, texturecoords,
                shape.quads_positions, shape.quads_normals,
                shape.quads_texturecoords, shape.positions, shape.normals,
                shape.texturecoords);
            if (!positions.empty())
                init_opengl_array_buffer(
                    vbos.positions_buffer, positions, false);
            if (!normals.empty())
                init_opengl_array_buffer(vbos.normals_buffer, normals, false);
            if (!texturecoords.empty())
                init_opengl_array_buffer(
                    vbos.texturecoords_buffer, texturecoords, false);
            if (!quads.empty()) {
                auto triangles = vector<vec3i>{};
                convert_quads_to_triangles(triangles, quads);
                init_opengl_elementbuffer(vbos.quads_buffer, triangles, false);
            }
        }
        state.shapes[shape_id] = vbos;
    }
}

void delete_drawgl_shape(drawgl_shape& glshape) {
    delete_opengl_array_buffer(glshape.positions_buffer);
    delete_opengl_array_buffer(glshape.normals_buffer);
    delete_opengl_array_buffer(glshape.texturecoords_buffer);
    delete_opengl_array_buffer(glshape.colors_buffer);
    delete_opengl_array_buffer(glshape.tangentspaces_buffer);
    delete_opengl_elementbuffer(glshape.points_buffer);
    delete_opengl_elementbuffer(glshape.lines_buffer);
    delete_opengl_elementbuffer(glshape.triangles_buffer);
    delete_opengl_elementbuffer(glshape.quads_buffer);
}

// delete state
void delete_drawgl_state(drawgl_state& state) {
    if (!state.program) return;
    delete_opengl_program(state.program);
    for (auto& texture : state.textures) delete_opengl_texture(texture);
    for (auto& shape : state.shapes) delete_drawgl_shape(shape);
    state.textures.clear();
    state.shapes.clear();
}

// draw with shading
void draw_widgets(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);

    begin_opengl_widgets_frame(win);
    if (begin_opengl_widgets_window(win, "yview")) {
        if (begin_tabbar_opengl_widget(win, "tabs")) {
            if (begin_tabitem_opengl_widget(win, "scene")) {
                draw_label_opengl_widget(
                    win, "scene", get_filename(app.filename));
                if (draw_button_opengl_widget(win, "load")) {
                    delete_drawgl_state(app.state);
                    load_scene_async(app);
                }
                draw_label_opengl_widget(win, "filename", app.filename);
                draw_label_opengl_widget(win, "status", app.status);
                end_tabitem_opengl_widget(win);
            }
            if (begin_tabitem_opengl_widget(win, "view")) {
                if (app.load_done) {
                    draw_combobox_opengl_widget(win, "camera",
                        app.draw_options.camera_id, app.scene.cameras, false);
                }
                draw_slider_opengl_widget(
                    win, "width", app.draw_options.image_width, 0, 4096);
                draw_slider_opengl_widget(
                    win, "height", app.draw_options.image_height, 0, 4096);
                draw_checkbox_opengl_widget(
                    win, "eyelight", app.draw_options.eyelight);
                continue_opengl_widget_line(win);
                draw_checkbox_opengl_widget(
                    win, "wireframe", app.draw_options.wireframe);
                continue_opengl_widget_line(win);
                draw_checkbox_opengl_widget(
                    win, "edges", app.draw_options.edges);
                if (app.time_range != zero2f) {
                    draw_slider_opengl_widget(win, "time", app.time,
                        app.time_range.x, app.time_range.y);
                    draw_textinput_opengl_widget(
                        win, "anim group", app.anim_group);
                    draw_checkbox_opengl_widget(win, "animate", app.animate);
                }
                draw_slider_opengl_widget(
                    win, "exposure", app.draw_options.exposure, -10, 10);
                draw_slider_opengl_widget(
                    win, "gamma", app.draw_options.gamma, 0.1f, 4);
                draw_checkbox_opengl_widget(
                    win, "double sided", app.draw_options.double_sided);
                draw_slider_opengl_widget(
                    win, "near", app.draw_options.near_plane, 0.01f, 1.0f);
                draw_slider_opengl_widget(
                    win, "far", app.draw_options.far_plane, 1000.0f, 10000.0f);
                draw_checkbox_opengl_widget(win, "fps", app.navigation_fps);
                if (draw_button_opengl_widget(win, "print cams")) {
                    for (auto& camera : app.scene.cameras) {
                        print_obj_camera(camera);
                    }
                }
                end_tabitem_opengl_widget(win);
            }
            if (app.load_done && begin_tabitem_opengl_widget(win, "navigate")) {
                draw_opengl_widgets_scene_tree(
                    win, "", app.scene, app.selection, app.update_list, 200);
                end_tabitem_opengl_widget(win);
            }
            if (app.load_done && begin_tabitem_opengl_widget(win, "inspect")) {
                draw_opengl_widgets_scene_inspector(
                    win, "", app.scene, app.selection, app.update_list, 200);
                end_tabitem_opengl_widget(win);
            }
            end_tabbar_opengl_widget(win);
        }
    }
    end_opengl_widgets_frame(win);
}

// draw with shading
void draw(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);

    clear_opengl_lframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.15f});
    set_opengl_viewport(get_opengl_framebuffer_size(win));
    if (app.load_done) {
        draw_glscene(app.state, app.scene, get_opengl_framebuffer_size(win),
            app.selection, app.draw_options);
    }
    draw_widgets(win);
    swap_opengl_buffers(win);
}

// update
void update(app_state& app) {
    // skip if not needed
    if (!app.load_done) return;

    // initialize gl state if needed
    if (!app.state.program) init_drawgl_state(app.state, app.scene);

    static auto last_time = 0.0f;
    for (auto& [type, index] : app.update_list) {
        if (type == typeid(yocto_texture)) {
            // TODO: update texture
            throw runtime_error("texture update not supported\n");
        }
        if (type == typeid(yocto_shape)) {
            // TODO: update shape
            throw runtime_error("shape update not supported\n");
        }
        if (type == typeid(yocto_scene_node) ||
            type == typeid(yocto_animation) || app.time != last_time) {
            update_transforms(app.scene, app.time, app.anim_group);
            last_time = app.time;
        }
    }
    app.update_list.clear();
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    delete_drawgl_state(app.state);
    app.filename = paths.front();
    load_scene_async(app);
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(
        win, {1280, 720}, "yview | " + get_filename(app.filename), &app, draw);

    // init widget
    init_opengl_widgets(win);

    // load textures and vbos
    update_transforms(app.scene, app.time);

    // loop
    auto mouse_pos = zero2f, last_pos = zero2f;
    auto last_time = get_time();
    while (!should_opengl_window_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_opengl_mouse_pos(win);
        auto mouse_left     = get_opengl_mouse_left(win);
        auto mouse_right    = get_opengl_mouse_right(win);
        auto alt_down       = get_opengl_alt_key(win);
        auto shift_down     = get_opengl_shift_key(win);
        auto widgets_active = get_opengl_widgets_active(win);

        // handle mouse and keyboard for navigation
        if (app.load_done && (mouse_left || mouse_right) && !alt_down &&
            !widgets_active) {
            auto dolly  = 0.0f;
            auto pan    = zero2f;
            auto rotate = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto& camera = app.scene.cameras.at(app.draw_options.camera_id);
            update_camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            app.update_list.push_back(
                {typeid(yocto_camera), app.draw_options.camera_id});
        }

        // animation
        if (app.load_done && app.animate) {
            auto time = (double)(get_time() - last_time) / 1000000000.0;
            app.time += min(1 / 60.0f, (float)time);
            if (app.time < app.time_range.x || app.time > app.time_range.y)
                app.time = app.time_range.x;
            update_transforms(app.scene, app.time);
            last_time = get_time();
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        process_opengl_events(
            win, !((mouse_left || mouse_right) || widgets_active));
    }

    // clear
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // initialize app
    app_state app{};
    auto      no_parallel = false;

    // parse command line
    auto parser = CLI::App{"views scenes inteactively"};
    parser.add_option("--camera", app.draw_options.camera_id, "Camera index.");
    parser.add_option("--hres,-R", app.draw_options.image_width,
        "Image horizontal resolution.");
    parser.add_option("--vres,-r", app.draw_options.image_height,
        "Image vertical resolution.");
    parser.add_flag("--eyelight,!--no-eyelight,-c", app.draw_options.eyelight,
        "Eyelight rendering.");
    parser.add_flag("--double-sided,!--no-double-sided,-D", app.double_sided,
        "Double-sided rendering.");
    parser.add_flag("--parallel,!--no-parallel", no_parallel,
        "Disable parallel execution.");
    parser.add_option("--output-image,-o", app.imfilename, "Image filename");
    parser.add_option("scene", app.filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix parallel code
    if (no_parallel) {
        app.load_options.run_serially = true;
    }

    // load
    load_scene_async(app);

    // run ui
    run_ui(app);

    // done
    return 0;
}
