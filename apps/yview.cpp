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

#include "../yocto/yocto_image.h"
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

struct drawgl_shape {
    opengl_array_buffer  positions_buffer     = {};
    opengl_array_buffer  normals_buffer       = {};
    opengl_array_buffer  texcoords_buffer     = {};
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

void init_drawgl_lights(drawgl_lights& lights, const yocto_scene& scene) {
    lights = {};
    for (auto& instance : scene.instances) {
        if (instance.shape < 0) continue;
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[instance.material];
        if (material.emission != zero3f) continue;
        if (lights.positions.size() >= 16) break;
        auto bbox = compute_bounds(shape);
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
}

// Draw options
struct draw_scene_params {
    int   camera_id        = 0;
    int   image_width      = 1280;
    int   image_height     = 720;
    bool  wireframe        = false;
    bool  edges            = false;
    float edge_offset      = 0.01f;
    bool  eyelight         = false;
    float exposure         = 0;
    float gamma            = 2.2f;
    vec3f ambient          = {0, 0, 0};
    bool  double_sided     = true;
    bool  non_rigid_frames = true;
    float near_plane       = 0.01f;
    float far_plane        = 10000.0f;
};

// Equality operators
inline bool operator==(const draw_scene_params& a, const draw_scene_params& b) {
    return memcmp(&a, &b, sizeof(a)) == 0;
}
inline bool operator!=(const draw_scene_params& a, const draw_scene_params& b) {
    return memcmp(&a, &b, sizeof(a)) != 0;
}

// Application task
enum struct app_task_type {
    none,
    load_scene,
    load_element,
    apply_edit,
    save_image,
    save_scene,
    close_scene
};

struct app_task {
    app_task_type                  type;
    future<void>                   result;
    atomic<bool>                   stop;
    atomic<int>                    current;
    concurrent_queue<image_region> queue;
    app_edit                       edit;

    app_task(app_task_type type, const app_edit& edit = {})
        : type{type}, result{}, stop{false}, current{-1}, queue{}, edit{edit} {}
    ~app_task() {
        stop = true;
        if (result.valid()) {
            try {
                result.get();
            } catch (...) {
            }
        }
    }
};

// Application state
struct app_scene {
    // loading parameters
    string filename  = "scene.json";
    string imagename = "out.png";
    string outname   = "scene.json";
    string name      = "";

    // options
    load_scene_params load_params = {};
    save_scene_params save_params = {};
    draw_scene_params draw_params = {};

    // scene
    yocto_scene scene = {};

    // rendering state
    drawgl_state  state  = {};
    drawgl_lights lights = {};

    // view image
    bool   navigation_fps = false;
    float  time           = 0;
    string anim_group     = "";
    vec2f  time_range     = zero2f;
    bool   animate        = false;

    // tasks
    bool            load_done = false;
    deque<app_task> task_queue;
    app_selection   selection = {typeid(void), -1};
};

// Application state
struct app_state {
    // data
    deque<app_scene> scenes;
    int              selected = -1;
    deque<string>    errors;

    // default options
    load_scene_params load_params = {};
    save_scene_params save_params = {};
    draw_scene_params draw_params = {};
};

void add_new_scene(app_state& app, const string& filename) {
    auto& scn       = app.scenes.emplace_back();
    scn.filename    = filename;
    scn.imagename   = get_noextension(filename) + ".png";
    scn.outname     = get_noextension(filename) + ".edited.yaml";
    scn.name        = get_filename(scn.filename);
    scn.load_params = app.load_params;
    scn.save_params = app.save_params;
    scn.draw_params = app.draw_params;
    scn.task_queue.emplace_back(app_task_type::load_scene);
    app.selected = (int)app.scenes.size() - 1;
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

        uniform mat4 shape_xform;                    // shape transform
        uniform mat4 shape_xform_invtranspose;       // shape transform
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
            norm = (shape_xform_invtranspose * vec4(norm,0)).xyz;
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

        uniform mat4 shape_xform;              // shape transform
        uniform mat4 shape_xform_invtranspose; // shape transform

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
    const draw_scene_params& options) {
    auto& shape    = scene.shapes[instance.shape];
    auto& vbos     = state.shapes.at(instance.shape);
    auto& material = scene.materials[instance.material];

    set_opengl_uniform(state.program, "shape_xform", mat4f(instance.frame));
    set_opengl_uniform(state.program, "shape_xform_invtranspose",
        transpose(mat4f(inverse(instance.frame, options.non_rigid_frames))));
    set_opengl_uniform(state.program, "shape_normal_offset", 0.0f);
    set_opengl_uniform(
        state.program, "highlight", (highlighted) ? vec4f{1, 1, 0, 1} : zero4f);

    auto mtype = 2;
    if (material.gltf_textures) mtype = 3;
    set_opengl_uniform(state.program, "mat_type", mtype);
    set_opengl_uniform(state.program, "mat_ke", material.emission);
    set_opengl_uniform(state.program, "mat_kd", material.diffuse);
    set_opengl_uniform(state.program, "mat_ks", vec3f{material.metallic});
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
        material.metallic_texture >= 0
            ? state.textures.at(material.metallic_texture)
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
        state.program, "vert_texcoord", vbos.texcoords_buffer, zero2f);
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
    const vec4i& viewport, const app_selection& highlighted,
    const draw_scene_params& options) {
    auto& camera      = scene.cameras.at(options.camera_id);
    auto  camera_view = mat4f(inverse(camera.frame));
    auto  camera_proj = make_perspective_mat(
        camera_fovx(camera) * (float)viewport.w / (float)viewport.z,
        (float)viewport.z / (float)viewport.w, options.near_plane,
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
            auto& material = scene.materials[instance.material];
            if (material.emission == zero3f) continue;
            if (lights_pos.size() >= 16) break;
            auto bbox = compute_bounds(shape);
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
        auto highlight = highlighted.type == typeid(yocto_instance) &&
                         highlighted.index == instance_id;
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
                true, true, true);
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
            if (!shape.texcoords.empty())
                init_opengl_array_buffer(
                    vbos.texcoords_buffer, shape.texcoords, false);
            if (!shape.colors.empty())
                init_opengl_array_buffer(
                    vbos.colors_buffer, shape.colors, false);
            if (!shape.tangents.empty())
                init_opengl_array_buffer(
                    vbos.tangentspaces_buffer, shape.tangents, false);
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
                quads_to_triangles(triangles, shape.quads);
                init_opengl_elementbuffer(vbos.quads_buffer, triangles, false);
            }
        } else {
            auto quads     = vector<vec4i>{};
            auto positions = vector<vec3f>{};
            auto normals   = vector<vec3f>{};
            auto texcoords = vector<vec2f>{};
            split_facevarying(quads, positions, normals, texcoords,
                shape.quads_positions, shape.quads_normals,
                shape.quads_texcoords, shape.positions, shape.normals,
                shape.texcoords);
            if (!positions.empty())
                init_opengl_array_buffer(
                    vbos.positions_buffer, positions, false);
            if (!normals.empty())
                init_opengl_array_buffer(vbos.normals_buffer, normals, false);
            if (!texcoords.empty())
                init_opengl_array_buffer(
                    vbos.texcoords_buffer, texcoords, false);
            if (!quads.empty()) {
                auto triangles = vector<vec3i>{};
                quads_to_triangles(triangles, quads);
                init_opengl_elementbuffer(vbos.quads_buffer, triangles, false);
            }
        }
        state.shapes[shape_id] = vbos;
    }
}

void delete_drawgl_shape(drawgl_shape& glshape) {
    delete_opengl_array_buffer(glshape.positions_buffer);
    delete_opengl_array_buffer(glshape.normals_buffer);
    delete_opengl_array_buffer(glshape.texcoords_buffer);
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
void draw_opengl_widgets(const opengl_window& win) {
    static string load_path = "", save_path = "", error_message = "";
    auto&         app = *(app_state*)get_opengl_user_pointer(win);
    if (!begin_opengl_widgets_window(win, "yview")) return;
    if (!app.errors.empty() && error_message.empty()) {
        error_message = app.errors.front();
        app.errors.pop_front();
        open_modal_opengl_widget(win, "error");
    }
    if (!draw_modal_message_opengl_window(win, "error", error_message)) {
        error_message = "";
    }
    if (draw_modal_fileialog_opengl_widgets(
            win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
        add_new_scene(app, load_path);
    }
    if (draw_modal_fileialog_opengl_widgets(win, "save", save_path, true,
            get_dirname(save_path), get_filename(save_path),
            "*.yaml;*.obj;*.pbrt")) {
        app.scenes[app.selected].outname = save_path;
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::save_scene);
        save_path = "";
    }
    if (draw_modal_fileialog_opengl_widgets(win, "save image", save_path, true,
            get_dirname(save_path), get_filename(save_path),
            "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
        app.scenes[app.selected].imagename = save_path;
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::save_image);
        save_path = "";
    }
    if (draw_button_opengl_widget(win, "load")) {
        open_modal_opengl_widget(win, "load");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save",
            app.selected >= 0 && app.scenes[app.selected].task_queue.empty())) {
        save_path = app.scenes[app.selected].outname;
        open_modal_opengl_widget(win, "save");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save image", app.selected >= 0)) {
        save_path = app.scenes[app.selected].imagename;
        open_modal_opengl_widget(win, "save image");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "close", app.selected >= 0)) {
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::close_scene);
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "quit")) {
        set_close_opengl_window(win, true);
    }
    if (app.scenes.empty()) return;
    draw_combobox_opengl_widget(
        win, "scene", app.selected, (int)app.scenes.size(),
        [&app](int idx) { return app.scenes[idx].name.c_str(); }, false);
    auto& scn = app.scenes[app.selected];
    if (begin_header_opengl_widget(win, "trace")) {
        auto cam_names = vector<string>();
        for (auto& camera : scn.scene.cameras) cam_names.push_back(camera.uri);
        auto draw_params = scn.draw_params;
        if (scn.load_done) {
            if (draw_combobox_opengl_widget(
                    win, "camera", draw_params.camera_id, cam_names)) {
            }
        }
        draw_slider_opengl_widget(
            win, "width", draw_params.image_width, 0, 4096);
        draw_slider_opengl_widget(
            win, "height", draw_params.image_height, 0, 4096);
        draw_checkbox_opengl_widget(win, "eyelight", draw_params.eyelight);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "wireframe", draw_params.wireframe);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "edges", draw_params.edges);
        if (scn.time_range != zero2f) {
            draw_slider_opengl_widget(
                win, "time", scn.time, scn.time_range.x, scn.time_range.y);
            draw_textinput_opengl_widget(win, "anim group", scn.anim_group);
            draw_checkbox_opengl_widget(win, "animate", scn.animate);
        }
        draw_slider_opengl_widget(
            win, "exposure", draw_params.exposure, -10, 10);
        draw_slider_opengl_widget(win, "gamma", draw_params.gamma, 0.1f, 4);
        draw_checkbox_opengl_widget(
            win, "double sided", draw_params.double_sided);
        draw_slider_opengl_widget(
            win, "near", draw_params.near_plane, 0.01f, 1.0f);
        draw_slider_opengl_widget(
            win, "far", draw_params.far_plane, 1000.0f, 10000.0f);

        if (draw_params != scn.draw_params) {
            scn.task_queue.emplace_back(app_task_type::apply_edit,
                app_edit{typeid(draw_scene_params), -1, draw_params, false});
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "inspect")) {
        draw_label_opengl_widget(win, "scene", get_filename(scn.filename));
        draw_label_opengl_widget(win, "filename", scn.filename);
        draw_label_opengl_widget(win, "outname", scn.outname);
        draw_label_opengl_widget(win, "imagename", scn.imagename);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "fps", scn.navigation_fps);
        if (draw_button_opengl_widget(win, "print cams")) {
            for (auto& camera : scn.scene.cameras) {
                print_obj_camera(camera);
            }
        }
        continue_opengl_widget_line(win);
        if (draw_button_opengl_widget(win, "print stats")) {
            print_info(format_stats(scn.scene));
        }
        end_header_opengl_widget(win);
    }
    if (scn.load_done && begin_header_opengl_widget(win, "scene tree")) {
        draw_opengl_widgets_scene_tree(win, "", scn.scene, scn.selection, 200);
        end_header_opengl_widget(win);
    }
    if (scn.load_done && begin_header_opengl_widget(win, "scene object")) {
        auto edit = app_edit{};
        if (draw_opengl_widgets_scene_inspector(
                win, "", scn.scene, scn.selection, edit, 200)) {
            scn.task_queue.emplace_back(app_task_type::apply_edit, edit);
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "log")) {
        draw_log_opengl_widget(win);
        end_header_opengl_widget(win);
    }
}

// draw with shading
void draw(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);

    clear_opengl_framebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    set_opengl_viewport(get_opengl_framebuffer_viewport(win));
    if (!app.scenes.empty() && app.selected >= 0 &&
        app.scenes[app.selected].load_done) {
        auto& scn = app.scenes[app.selected];
        draw_glscene(scn.state, scn.scene, get_opengl_framebuffer_viewport(win),
            scn.selection, scn.draw_params);
    }
    begin_opengl_widgets_frame(win);
    draw_opengl_widgets(win);
    end_opengl_widgets_frame(win);
    swap_opengl_buffers(win);
}

void apply_edit(const string& filename, yocto_scene& scene,
    drawgl_lights& lights, draw_scene_params& draw_params, float time,
    const string& anim_group, bool& reload_element, const app_edit& edit) {
    static auto last_time             = 0.0f;
    bool        updated_hierarchy     = false;
    auto& [type, index, data, reload] = edit;
    if (type == typeid(yocto_camera)) {
        scene.cameras[index] = any_cast<yocto_camera>(data);
    } else if (type == typeid(yocto_texture)) {
        scene.textures[index] = any_cast<yocto_texture>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_voltexture)) {
        scene.voltextures[index] = any_cast<yocto_voltexture>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_shape)) {
        scene.shapes[index] = any_cast<yocto_shape>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_subdiv)) {
        // TODO: this needs more fixing?
        scene.subdivs[index] = any_cast<yocto_subdiv>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_material)) {
        scene.materials[index] = any_cast<yocto_material>(data);
    } else if (type == typeid(yocto_instance)) {
        scene.instances[index] = any_cast<yocto_instance>(data);
    } else if (type == typeid(yocto_animation)) {
        scene.animations[index] = any_cast<yocto_animation>(data);
        updated_hierarchy       = true;
    } else if (type == typeid(yocto_scene_node)) {
        scene.nodes[index] = any_cast<yocto_scene_node>(data);
        updated_hierarchy  = true;
    } else if (type == typeid(draw_scene_params)) {
        draw_params = any_cast<draw_scene_params>(data);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }
    if (updated_hierarchy || time != last_time) {
        update_transforms(scene, time, anim_group);
        last_time = time;
    }
}

// reload an element
void load_element(
    const string& filename, yocto_scene& scene, const app_edit& edit) {
    auto& [type, index, data, reload] = edit;

    if (type == typeid(yocto_texture)) {
        auto& texture = scene.textures[index];
        load_image(get_dirname(filename) + texture.uri, texture.hdr_image,
            texture.ldr_image);
    } else if (type == typeid(yocto_voltexture)) {
        auto& texture = scene.voltextures[index];
        load_volume(get_dirname(filename) + texture.uri, texture.volume_data);
    } else if (type == typeid(yocto_shape)) {
        auto& shape = scene.shapes[index];
        load_shape(get_dirname(filename) + shape.uri, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.quads_positions,
            shape.quads_normals, shape.quads_texcoords, shape.positions,
            shape.normals, shape.texcoords, shape.colors, shape.radius, false);
    } else if (type == typeid(yocto_subdiv)) {
        // TODO: this needs more fixing?
        auto& subdiv = scene.subdivs[index];
        load_shape(get_dirname(filename) + subdiv.uri, subdiv.points,
            subdiv.lines, subdiv.triangles, subdiv.quads,
            subdiv.quads_positions, subdiv.quads_normals,
            subdiv.quads_texcoords, subdiv.positions, subdiv.normals,
            subdiv.texcoords, subdiv.colors, subdiv.radius,
            subdiv.preserve_facevarying);
        tesselate_subdiv(scene, scene.subdivs[index]);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }
}

// update
void update(app_state& app) {
    // close if needed
    while (!app.scenes.empty()) {
        auto pos = -1;
        for (auto idx = 0; idx < app.scenes.size(); idx++) {
            for (auto& task : app.scenes[idx].task_queue) {
                if (task.type == app_task_type::close_scene) pos = idx;
            }
        }
        if (pos < 0) break;
        app.scenes.erase(app.scenes.begin() + pos);
        app.selected = app.scenes.empty() ? -1 : 0;
    }

    // apply synchronous edit
    for (auto& scn : app.scenes) {
        while (!scn.task_queue.empty()) {
            auto& task = scn.task_queue.front();
            if (task.type != app_task_type::apply_edit) break;
            log_info("start editing " + scn.filename);
            try {
                auto reload_element = false;
                apply_edit(scn.filename, scn.scene, scn.lights, scn.draw_params,
                    scn.time, scn.anim_group, reload_element, task.edit);
                log_info("done editing " + scn.filename);
                if (reload_element) {
                    scn.load_done = false;
                    scn.task_queue.emplace_back(
                        app_task_type::load_element, task.edit);
                }
            } catch (std::exception& e) {
                log_error(e.what());
                app.errors.push_back("cannot edit " + scn.filename);
            }
            scn.task_queue.pop_front();
        }
    }

    // grab result of finished tasks
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (!task.result.valid() || !task.queue.empty() ||
            task.result.wait_for(std::chrono::nanoseconds(10)) !=
                std::future_status::ready)
            continue;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::close_scene: break;
            case app_task_type::load_scene: {
                try {
                    task.result.get();
                    scn.load_done = true;
                    init_drawgl_state(scn.state, scn.scene);
                    log_info("done loading " + scn.filename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = get_filename(scn.filename) + " [error]";
                    app.errors.push_back("cannot load " + scn.filename);
                }
            } break;
            case app_task_type::load_element: {
                try {
                    task.result.get();
                    scn.load_done = true;
                    log_info("done loading element from " + scn.filename);
                    if (task.edit.type == typeid(yocto_texture)) {
                        // not supported yet
                        log_error("texture refresh is not supported yet");
                    } else if (task.edit.type == typeid(yocto_voltexture)) {
                    } else if (task.edit.type == typeid(yocto_shape)) {
                        // not supported yet
                        log_error("shape refresh is not supported yet");
                    } else if (task.edit.type == typeid(yocto_subdiv)) {
                        // not supported yet
                        log_error("shape refresh is not supported yet");
                    } else {
                        throw runtime_error("unsupported type");
                    }
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = get_filename(scn.filename) + " [error]";
                    app.errors.push_back(
                        "cannot load element from " + scn.filename);
                }
            } break;
            case app_task_type::save_image: {
                try {
                    task.result.get();
                    log_info("done saving " + scn.imagename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot save " + scn.imagename);
                }
            } break;
            case app_task_type::save_scene: {
                try {
                    task.result.get();
                    log_info("done saving " + scn.outname);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot save " + scn.outname);
                }
            } break;
            case app_task_type::apply_edit: break;
        }
        scn.task_queue.pop_front();
    }
    // schedule tasks not running
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (task.result.valid()) continue;
        task.stop = false;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::close_scene: break;
            case app_task_type::load_scene: {
                log_info("start loading " + scn.filename);
                scn.load_done = false;
                task.result   = async([&scn]() {
                    load_scene(scn.filename, scn.scene, scn.load_params);
                    tesselate_subdivs(scn.scene);
                    init_drawgl_lights(scn.lights, scn.scene);
                    if (scn.lights.empty() && !scn.draw_params.eyelight) {
                        print_info(
                            "no lights presents, switching to eyelight shader");
                        scn.draw_params.eyelight = true;
                    }
                    scn.time_range = compute_animation_range(scn.scene);
                    scn.time       = scn.time_range.x;
                });
            } break;
            case app_task_type::load_element: {
                log_info("start loading element for " + scn.filename);
                scn.load_done = false;
                task.result   = async([&scn, &task]() {
                    load_element(scn.filename, scn.scene, task.edit);
                });
            } break;
            case app_task_type::save_image: {
                log_info("start saving " + scn.imagename);
                task.result = async(
                    []() { throw runtime_error("not implemnted yet"); });
            } break;
            case app_task_type::save_scene: {
                log_info("start saving " + scn.outname);
                task.result = async([&scn]() {
                    save_scene(scn.outname, scn.scene, scn.save_params);
                });
            } break;
            case app_task_type::apply_edit: break;
        }
    }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    for (auto& path : paths) add_new_scene(app, path);
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280 + 320, 720}, "yview", &app, draw);

    // init widget
    init_opengl_widgets(win);

    // setup logging
    set_log_callback(
        [&win](const string& msg) { add_log_opengl_widget(win, msg.c_str()); });

    // loop
    auto mouse_pos = zero2f, last_pos = zero2f;
    auto last_time = std::chrono::high_resolution_clock::now();
    while (!should_opengl_window_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_opengl_mouse_pos(win);
        auto mouse_left     = get_opengl_mouse_left(win);
        auto mouse_right    = get_opengl_mouse_right(win);
        auto alt_down       = get_opengl_alt_key(win);
        auto shift_down     = get_opengl_shift_key(win);
        auto widgets_active = get_opengl_widgets_active(win);

        // update trasforms
        if (app.selected >= 0) {
            auto& scn = app.scenes[app.selected];
            update_transforms(scn.scene, scn.time);
        }

        // handle mouse and keyboard for navigation
        if (app.selected >= 0 && (mouse_left || mouse_right) && !alt_down &&
            !widgets_active) {
            auto& scn        = app.scenes[app.selected];
            auto& old_camera = scn.scene.cameras.at(scn.draw_params.camera_id);
            auto  camera     = scn.scene.cameras.at(scn.draw_params.camera_id);
            auto  dolly      = 0.0f;
            auto  pan        = zero2f;
            auto  rotate     = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            update_camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            if (camera.frame != old_camera.frame ||
                camera.focus_distance != old_camera.focus_distance) {
                scn.task_queue.emplace_back(app_task_type::apply_edit,
                    app_edit{typeid(yocto_camera), scn.draw_params.camera_id,
                        camera, false});
            }
        }

        // animation
        if (app.selected >= 0 && app.scenes[app.selected].animate) {
            auto& scn     = app.scenes[app.selected];
            auto  now     = std::chrono::high_resolution_clock::now();
            auto  elapsed = now - last_time;
            auto  time    = (double)(elapsed.count()) / 1000000000.0;
            scn.time += min(1 / 60.0f, (float)time);
            if (scn.time < scn.time_range.x || scn.time > scn.time_range.y)
                scn.time = scn.time_range.x;
            update_transforms(scn.scene, scn.time);
            last_time = now;
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        process_opengl_events(win);
    }

    // clear
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // initialize app
    app_state app{};
    auto      filenames   = vector<string>{};
    auto      no_parallel = false;

    // parse command line
    auto parser = CLI::App{"views scenes inteactively"};
    parser.add_option("--camera", app.draw_params.camera_id, "Camera index.");
    parser.add_option("--hres,-R", app.draw_params.image_width,
        "Image horizontal resolution.");
    parser.add_option("--vres,-r", app.draw_params.image_height,
        "Image vertical resolution.");
    parser.add_flag("--eyelight,!--no-eyelight,-c", app.draw_params.eyelight,
        "Eyelight rendering.");
    parser.add_flag("--parallel,!--no-parallel", no_parallel,
        "Disable parallel execution.");
    parser.add_option("scenes", filenames, "Scene filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix parallel code
    if (no_parallel) {
        app.load_params.run_serially = true;
        app.save_params.run_serially = true;
    }

    // loading images
    for (auto filename : filenames) add_new_scene(app, filename);
    app.selected = app.scenes.empty() ? -1 : 0;

    // run ui
    run_ui(app);

    // done
    return 0;
}
