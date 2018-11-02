//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"
#include "yglutils.h"
#include "ysceneui.h"

struct glshape {
    glarraybuffer gl_pos = {}, gl_norm = {}, gl_texcoord = {}, gl_color = {},
                  gl_tangsp   = {};
    glelementbuffer gl_points = {}, gl_lines = {}, gl_triangles = {},
                    gl_quads                      = {};
    vector<glelementbuffer> gl_split_quads        = {};
    int                     num_facevarying_quads = 0;
};

struct draw_glstate {
    glprogram         prog     = {};
    vector<glshape>   shapes   = {};
    vector<glshape>   surfaces = {};
    vector<gltexture> textures = {};
};

// Application state
struct app_state {
    // scene
    yocto_scene scene = {};

    // parameters
    string filename    = "scene.json";  // scene name
    string imfilename  = "out.png";     // output image
    string outfilename = "scene.json";  // save scene name
    int    camid       = 0;             // camera id
    int    resolution  = 512;           // image resolution
    bool   wireframe   = false;         // wireframe drawing
    bool   edges       = false;         // draw edges
    float  edge_offset = 0.01f;         // offset for edges
    bool   eyelight    = false;         // camera light mode
    float  exposure    = 0;             // exposure
    float  gamma       = 2.2f;          // gamma
    vec3f  ambient     = {0, 0, 0};     // ambient lighting
    float  near_plane  = 0.01f;         // near plane
    float  far_plane   = 10000.0f;      // far plane

    draw_glstate state = {};

    bool                       widgets_open   = false;
    bool                       navigation_fps = false;
    tuple<string, int>         selection      = {"", -1};
    vector<tuple<string, int>> update_list;
    float                      time       = 0;
    string                     anim_group = "";
    vec2f                      time_range = zero2f;
    bool                       animate    = false;
};

void draw_glscene(draw_glstate& state, const yocto_scene& scene,
    const yocto_camera& camera, const vec2i& viewport_size,
    const tuple<string, int>& highlighted, bool eyelight, bool wireframe,
    bool edges, float exposure, float gamma, float near_plane, float far_plane);

// draw with shading
void draw(const glwindow& win) {
    auto& app              = *(app_state*)get_user_pointer(win);
    auto  framebuffer_size = get_glframebuffer_size(win);
    app.resolution         = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app.update_list) {
        if (get<0>(sel) == "texture") {
            // TODO: update texture
            printf("texture update not supported\n");
        }
        if (get<0>(sel) == "subdiv") {
            // TODO: update subdiv
            printf("subdiv update not supported\n");
        }
        if (get<0>(sel) == "shape") {
            // TODO: update shape
            printf("shape update not supported\n");
        }
        if (get<0>(sel) == "node" || get<0>(sel) == "animation" ||
            app.time != last_time) {
            update_transforms(app.scene, app.time, app.anim_group);
            last_time = app.time;
        }
    }
    app.update_list.clear();

    auto& camera = app.scene.cameras.at(app.camid);
    clear_glframebuffer(vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    draw_glscene(app.state, app.scene, camera, framebuffer_size, app.selection,
        app.eyelight, app.wireframe, app.edges, app.exposure, app.gamma,
        app.near_plane, app.far_plane);

    begin_glwidgets_frame(win);
    if (begin_glwidgets_window(win, "yview")) {
        if (begin_header_glwidget(win, "scene")) {
            draw_label_glwidgets(win, "scene", "%s", app.filename.c_str());
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "view")) {
            draw_combobox_glwidget(
                win, "camera", app.camid, app.scene.cameras, false);
            draw_slider_glwidget(win, "resolution", app.resolution, 256, 4096);
            draw_checkbox_glwidget(win, "eyelight", app.eyelight);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "wireframe", app.wireframe);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "edges", app.edges);
            if (app.time_range != zero2f) {
                draw_slider_glwidget(
                    win, "time", app.time, app.time_range.x, app.time_range.y);
                draw_textinput_glwidget(win, "anim group", app.anim_group);
                draw_checkbox_glwidget(win, "animate", app.animate);
            }
            draw_slider_glwidget(win, "exposure", app.exposure, -10, 10);
            draw_slider_glwidget(win, "gamma", app.gamma, 0.1f, 4);
            draw_slider_glwidget(win, "near", app.near_plane, 0.01f, 1.0f);
            draw_slider_glwidget(win, "far", app.far_plane, 1000.0f, 10000.0f);
            draw_checkbox_glwidget(win, "fps", app.navigation_fps);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "navigate")) {
            draw_glwidgets_scene_tree(
                win, "", app.scene, app.selection, app.update_list, 200);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "inspect")) {
            draw_glwidgets_scene_inspector(
                win, "", app.scene, app.selection, app.update_list, 200);
            end_header_glwidget(win);
        }
    }
    end_glwidgets_frame(win);

    swap_glbuffers(win);
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
void draw_glinstance(draw_glstate& state, const yocto_scene& scene,
    const yocto_instance& instance, bool highlighted, bool eyelight, bool edges) {
    if (instance.shape >= 0) {
        auto& shape    = scene.shapes[instance.shape];
        auto& vbos     = state.shapes.at(instance.shape);
        auto& material = scene.materials[shape.material];

        auto xform = frame_to_mat(instance.frame);

        set_gluniform(state.prog, "shape_xform", xform);
        set_gluniform(state.prog, "shape_normal_offset", 0.0f);

        auto mtype = 1;
        if (material.base_metallic) mtype = 2;
        if (material.gltf_textures) mtype = (material.base_metallic) ? 2 : 3;
        set_gluniform(state.prog, "mat_type", mtype);
        set_gluniform(state.prog, "mat_ke", material.emission);
        set_gluniform(state.prog, "mat_kd", material.diffuse);
        set_gluniform(state.prog, "mat_ks", material.specular);
        set_gluniform(state.prog, "mat_rs", material.roughness);
        set_gluniform(state.prog, "mat_op", material.opacity);
        set_gluniform(state.prog, "mat_double_sided", (int)material.double_sided);
        set_gluniform_texture(state.prog, "mat_ke_txt", "mat_ke_txt_on",
            material.emission_texture >= 0 ?
                state.textures.at(material.emission_texture) :
                gltexture{},
            0);
        set_gluniform_texture(state.prog, "mat_kd_txt", "mat_kd_txt_on",
            material.diffuse_texture >= 0 ?
                state.textures.at(material.diffuse_texture) :
                gltexture{},
            1);
        set_gluniform_texture(state.prog, "mat_ks_txt", "mat_ks_txt_on",
            material.specular_texture >= 0 ?
                state.textures.at(material.specular_texture) :
                gltexture{},
            2);
        set_gluniform_texture(state.prog, "mat_rs_txt", "mat_rs_txt_on",
            material.roughness_texture >= 0 ?
                state.textures.at(material.roughness_texture) :
                gltexture{},
            3);
        set_gluniform_texture(state.prog, "mat_op_txt", "mat_op_txt_on",
            material.opacity_texture >= 0 ?
                state.textures.at(material.opacity_texture) :
                gltexture{},
            4);
        set_gluniform_texture(state.prog, "mat_norm_txt", "mat_norm_txt_on",
            material.normal_texture >= 0 ?
                state.textures.at(material.normal_texture) :
                gltexture{},
            5);

        set_gluniform(state.prog, "elem_faceted", (int)shape.normals.empty());
        set_glvertexattrib(state.prog, "vert_pos", vbos.gl_pos, zero3f);
        set_glvertexattrib(state.prog, "vert_norm", vbos.gl_norm, zero3f);
        set_glvertexattrib(state.prog, "vert_texcoord", vbos.gl_texcoord, zero2f);
        set_glvertexattrib(
            state.prog, "vert_color", vbos.gl_color, vec4f{1, 1, 1, 1});
        set_glvertexattrib(
            state.prog, "vert_tangsp", vbos.gl_tangsp, vec4f{0, 0, 1, 1});

        if (vbos.gl_points) {
            set_gluniform(state.prog, "elem_type", 1);
            draw_glpoints(vbos.gl_points, vbos.gl_points.num);
        }
        if (vbos.gl_lines) {
            set_gluniform(state.prog, "elem_type", 2);
            draw_gllines(vbos.gl_lines, vbos.gl_lines.num);
        }
        if (vbos.gl_triangles) {
            set_gluniform(state.prog, "elem_type", 3);
            draw_gltriangles(vbos.gl_triangles, vbos.gl_triangles.num);
        }
        if (vbos.gl_quads) {
            set_gluniform(state.prog, "elem_type", 3);
            draw_gltriangles(vbos.gl_quads, vbos.gl_quads.num);
        }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_glerror();
        set_gluniform(state.prog, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state.prog, "ke"), 0, 0, 0);
        set_gluniform(state.prog, "op"), material.op);
        set_gluniform(state.prog, "shp_normal_offset"), 0.01f);
        check_glerror();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_glerror();
    }
#endif
        if (edges) printf("edges are momentarily disabled\n");

        // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
    } else if (instance.surface >= 0) {
        auto& surface = scene.surfaces[instance.surface];
        auto& vbos    = state.surfaces.at(instance.surface);
        for (auto group_id = 0; group_id < vbos.gl_split_quads.size();
             group_id++) {
            auto& material = scene.materials[surface.materials.at(group_id)];

            auto xform = frame_to_mat(instance.frame);

            set_gluniform(state.prog, "shape_xform", xform);
            set_gluniform(state.prog, "shape_normal_offset", 0.0f);

            auto mtype = 1;
            if (material.base_metallic) mtype = 2;
            if (material.gltf_textures)
                mtype = (material.base_metallic) ? 2 : 3;
            set_gluniform(state.prog, "mat_type", mtype);
            set_gluniform(state.prog, "mat_ke", material.emission);
            set_gluniform(state.prog, "mat_kd", material.diffuse);
            set_gluniform(state.prog, "mat_ks", material.specular);
            set_gluniform(state.prog, "mat_rs", material.roughness);
            set_gluniform(state.prog, "mat_op", material.opacity);
            set_gluniform(
                state.prog, "mat_double_sided", (int)material.double_sided);
            set_gluniform_texture(state.prog, "mat_ke_txt", "mat_ke_txt_on",
                material.emission_texture >= 0 ?
                    state.textures.at(material.emission_texture) :
                    gltexture{},
                0);
            set_gluniform_texture(state.prog, "mat_kd_txt", "mat_kd_txt_on",
                material.diffuse_texture >= 0 ?
                    state.textures.at(material.diffuse_texture) :
                    gltexture{},
                1);
            set_gluniform_texture(state.prog, "mat_ks_txt", "mat_ks_txt_on",
                material.specular_texture >= 0 ?
                    state.textures.at(material.specular_texture) :
                    gltexture{},
                2);
            set_gluniform_texture(state.prog, "mat_rs_txt", "mat_rs_txt_on",
                material.roughness_texture >= 0 ?
                    state.textures.at(material.roughness_texture) :
                    gltexture{},
                3);
            set_gluniform_texture(state.prog, "mat_op_txt", "mat_op_txt_on",
                material.opacity_texture >= 0 ?
                    state.textures.at(material.opacity_texture) :
                    gltexture{},
                4);
            set_gluniform_texture(state.prog, "mat_norm_txt", "mat_norm_txt_on",
                material.normal_texture >= 0 ?
                    state.textures.at(material.normal_texture) :
                    gltexture{},
                5);

            set_gluniform(
                state.prog, "elem_faceted", (int)surface.normals.empty());
            set_glvertexattrib(state.prog, "vert_pos", vbos.gl_pos, zero3f);
            set_glvertexattrib(state.prog, "vert_norm", vbos.gl_norm, zero3f);
            set_glvertexattrib(
                state.prog, "vert_texcoord", vbos.gl_texcoord, zero2f);
            set_glvertexattrib(
                state.prog, "vert_color", vbos.gl_color, vec4f{1, 1, 1, 1});
            set_glvertexattrib(
                state.prog, "vert_tangsp", vbos.gl_tangsp, vec4f{0, 0, 1, 1});

            if (vbos.gl_split_quads[group_id]) {
                set_gluniform(state.prog, "elem_type", 3);
                draw_gltriangles(vbos.gl_split_quads[group_id],
                    vbos.gl_split_quads[group_id].num);
            }
        }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_glerror();
        set_gluniform(state.prog, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state.prog, "ke"), 0, 0, 0);
        set_gluniform(state.prog, "op"), material.op);
        set_gluniform(state.prog, "shp_normal_offset"), 0.01f);
        check_glerror();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_glerror();
    }
#endif
        if (edges) printf("edges are momentarily disabled\n");

        // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
    }
}

// Display a scene
void draw_glscene(draw_glstate& state, const yocto_scene& scene,
    const yocto_camera& camera, const vec2i& viewport_size,
    const tuple<string, int>& highlighted, bool eyelight, bool wireframe,
    bool edges, float exposure, float gamma, float near_plane, float far_plane) {
    set_glviewport(viewport_size);

    auto camera_view = frame_to_mat(inverse(camera.frame));
    // auto camera_proj =
    //         perspective_mat(evaluate_camera_fovy(camera),
    //             (float)viewport_size.x / (float)viewport_size.y, near_plane);
    auto camera_proj = perspective_mat(evaluate_camera_fovy(camera),
        (float)viewport_size.x / (float)viewport_size.y, near_plane, far_plane);

    bind_glprogram(state.prog);
    set_gluniform(state.prog, "cam_pos", camera.frame.o);
    set_gluniform(state.prog, "cam_xform_inv", camera_view);
    set_gluniform(state.prog, "cam_proj", camera_proj);
    set_gluniform(state.prog, "eyelight", (int)eyelight);
    set_gluniform(state.prog, "exposure", exposure);
    set_gluniform(state.prog, "gamma", gamma);

    if (!eyelight) {
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
                    area += quad_area(shape.positions[q.x], shape.positions[q.y],
                        shape.positions[q.z], shape.positions[q.w]);
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
        if (lights_pos.empty()) eyelight = false;
        set_gluniform(state.prog, "lamb", zero3f);
        set_gluniform(state.prog, "lnum", (int)lights_pos.size());
        for (auto i = 0; i < lights_pos.size(); i++) {
            auto is = to_string(i);
            set_gluniform(
                state.prog, ("lpos[" + is + "]").c_str(), lights_pos[i]);
            set_gluniform(state.prog, ("lke[" + is + "]").c_str(), lights_ke[i]);
            set_gluniform(
                state.prog, ("ltype[" + is + "]").c_str(), (int)lights_type[i]);
        }
    }

    if (wireframe) set_glwireframe(true);
    for (auto instance_id = 0; instance_id < scene.instances.size();
         instance_id++) {
        auto& instance = scene.instances[instance_id];
        // auto& shape     = scene.shapes[instance.shape];
        // auto& material  = scene.materials[shape.material];
        auto highlight = highlighted ==
                         tuple<string, int>{"instance", instance_id};
        draw_glinstance(state, scene, instance, highlight, eyelight, edges);
    }

    unbind_glprogram();
    if (wireframe) set_glwireframe(false);
}

draw_glstate init_draw_state(const glwindow& win) {
    auto& app   = *(app_state*)get_user_pointer(win);
    auto  state = draw_glstate();
    // load textures and vbos
    state.prog = make_glprogram(vertex, fragment);
    state.textures.resize(app.scene.textures.size());
    for (auto texture_id = 0; texture_id < app.scene.textures.size();
         texture_id++) {
        auto texture = app.scene.textures[texture_id];
        if (!texture.hdr_image.pixels.empty()) {
            state.textures[texture_id] = make_gltexture(
                texture.hdr_image, true, true, true);
        } else if (!texture.ldr_image.pixels.empty()) {
            state.textures[texture_id] = make_gltexture(
                texture.ldr_image, !texture.ldr_as_linear, true, true);
        } else {
            printf("bad texture");
        }
    }
    state.shapes.resize(app.scene.shapes.size());
    for (auto shape_id = 0; shape_id < app.scene.shapes.size(); shape_id++) {
        auto& shape = app.scene.shapes[shape_id];
        auto  vbos  = glshape();
        if (!shape.positions.empty())
            vbos.gl_pos = make_glarraybuffer(shape.positions, false);
        if (!shape.normals.empty())
            vbos.gl_norm = make_glarraybuffer(shape.normals, false);
        if (!shape.texturecoords.empty())
            vbos.gl_texcoord = make_glarraybuffer(shape.texturecoords, false);
        if (!shape.colors.empty())
            vbos.gl_color = make_glarraybuffer(shape.colors, false);
        if (!shape.tangentspaces.empty())
            vbos.gl_tangsp = make_glarraybuffer(shape.tangentspaces, false);
        if (!shape.points.empty())
            vbos.gl_points = make_glelementbuffer(shape.points, false);
        if (!shape.lines.empty())
            vbos.gl_lines = make_glelementbuffer(shape.lines, false);
        if (!shape.triangles.empty())
            vbos.gl_triangles = make_glelementbuffer(shape.triangles, false);
        if (!shape.quads.empty())
            vbos.gl_quads = make_glelementbuffer(
                convert_quads_to_triangles(shape.quads), false);
        state.shapes[shape_id] = vbos;
    }
    state.surfaces.resize(app.scene.surfaces.size());
    for (auto surface_id = 0; surface_id < app.scene.surfaces.size();
         surface_id++) {
        auto& surface       = app.scene.surfaces[surface_id];
        auto  vbos          = glshape();
        auto  quads         = vector<vec4i>();
        auto  positions     = vector<vec3f>();
        auto  normals       = vector<vec3f>();
        auto  texturecoords = vector<vec2f>();
        tie(quads, positions, normals, texturecoords) = convert_face_varying(
            surface.quads_positions, surface.quads_normals,
            surface.quads_texturecoords, surface.positions, surface.normals,
            surface.texturecoords);
        auto split_quads = vector<vector<vec4i>>();
        if (surface.materials.size() > 1 && !surface.quads_materials.empty()) {
            split_quads = ungroup_quads(quads, surface.quads_materials);
        } else {
            split_quads = {quads};
        }
        if (!positions.empty())
            vbos.gl_pos = make_glarraybuffer(positions, false);
        if (!normals.empty()) vbos.gl_norm = make_glarraybuffer(normals, false);
        if (!texturecoords.empty())
            vbos.gl_texcoord = make_glarraybuffer(texturecoords, false);
        vbos.gl_split_quads = {};
        for (auto& quads : split_quads) {
            if (!quads.empty()) vbos.gl_split_quads.push_back({});
            vbos.gl_split_quads.back() = make_glelementbuffer(
                convert_quads_to_triangles(quads), false);
        }
        vbos.num_facevarying_quads = (int)quads.size();
        state.surfaces[surface_id] = vbos;
    }
    return state;
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto& camera = app.scene.cameras.at(app.camid);
    auto width = clamp(evaluate_image_size(camera, app.resolution).x, 256, 1440),
         height = clamp(evaluate_image_size(camera, app.resolution).y, 256, 1440);
    auto win    = glwindow();
    init_glwindow(win, width, height, "yview", &app, draw);

    // init widget
    init_glwidgets(win);

    // load textures and vbos
    update_transforms(app.scene, app.time);

    // init gl data
    app.state = init_draw_state(win);

    // loop
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_glwindow_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_glmouse_pos(win);
        auto mouse_left     = get_glmouse_left(win);
        auto mouse_right    = get_glmouse_right(win);
        auto alt_down       = get_glalt_key(win);
        auto shift_down     = get_glshift_key(win);
        auto widgets_active = get_glwidgets_active(win);

        // handle mouse and keyboard for navigation
        if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
            auto dolly  = 0.0f;
            auto pan    = zero2f;
            auto rotate = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto& camera = app.scene.cameras.at(app.camid);
            camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            app.update_list.push_back({"camera", app.camid});
        }

        // animation
        if (app.animate) {
            app.time += 1 / 60.0f;
            if (app.time < app.time_range.x || app.time > app.time_range.y)
                app.time = app.time_range.x;
            update_transforms(app.scene, app.time);
        }

        // draw
        draw(win);

        // event hadling
        process_glevents(win, !((mouse_left || mouse_right) || widgets_active));
    }

    // clear
    delete_glwindow(win);
}

// Load INI file. The implementation does not handle escaping.
unordered_map<string, unordered_map<string, string>> load_ini(
    const string& filename) {
    auto f = fopen(filename.c_str(), "rt");
    if (!f) {
        log_error("cannot open {}", filename);
        return {};
    }
    auto ret       = unordered_map<string, unordered_map<string, string>>();
    auto cur_group = string();
    ret[""]        = {};

    char buf[4096];
    while (fgets(buf, 4096, f)) {
        auto line = string(buf);
        if (line.empty()) continue;
        if (line.front() == ';') continue;
        if (line.front() == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') {
                log_error("bad INI format");
                return {};
            }
            cur_group      = line.substr(1, line.length() - 2);
            ret[cur_group] = {};
        } else if (line.find('=') != line.npos) {
            auto var            = line.substr(0, line.find('='));
            auto val            = line.substr(line.find('=') + 1);
            ret[cur_group][var] = val;
        } else {
            log_error("bad INI format");
            return {};
        }
    }

    fclose(f);

    return ret;
}

int main(int argc, char* argv[]) {
    // initialize app
    auto app = app_state();

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "views scenes inteactively", "yview");
    app.camid      = parse_arg(parser, "--camera", 0, "Camera index.");
    app.resolution = parse_arg(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app.eyelight = parse_arg(
        parser, "--eyelight,-c", false, "Eyelight rendering.");
    auto double_sided = parse_arg(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto quiet = parse_arg(
        parser, "--quiet,-q", false, "Print only errors messages");
    auto highlight_filename = parse_arg(
        parser, "--highlights", ""s, "Highlight filename");
    app.imfilename = parse_arg(
        parser, "--output-image,-o", "out.png"s, "Image filename");
    app.filename = parse_arg(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    if (!load_scene(app.filename, app.scene))
        log_fatal("cannot load scene {}", app.filename);

    // tesselate
    if (!quiet) log_info("tesselating scene elements\n");
    tesselate_shapes_and_surfaces(app.scene);

    // add components
    if (!quiet) log_info("adding scene elements\n");
    if (double_sided) {
        for (auto& material : app.scene.materials) material.double_sided = true;
    }
    for (auto& err : validate_scene(app.scene)) log_error(err);

    // animation
    auto time_range = compute_animation_range(app.scene);
    app.time        = time_range.x;

    // run ui
    run_ui(app);

    // done
    return 0;
}
