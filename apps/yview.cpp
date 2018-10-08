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

struct draw_glshape_vbos {
    unsigned int gl_pos = 0, gl_norm = 0, gl_texcoord = 0, gl_color = 0,
                 gl_tangsp = 0, gl_points = 0, gl_lines = 0,
                 gl_triangles = 0;  // unmanaged data for OpenGL viewer
};

struct draw_glstate {
    unsigned int                                        gl_prog = 0;
    std::unordered_map<const shape*, draw_glshape_vbos> shp_vbos;
    std::unordered_map<const texture*, unsigned int>    txt_id;
};

// Application state
struct app_state {
    // scene
    scene* scn = nullptr;

    // parameters
    std::string filename    = "scene.json";  // scene name
    std::string imfilename  = "out.png";     // output image
    std::string outfilename = "scene.json";  // save scene name
    int         camid       = 0;             // camera id
    int         resolution  = 512;           // image resolution
    bool        wireframe   = false;         // wireframe drawing
    bool        edges       = false;         // draw edges
    float       edge_offset = 0.01f;         // offset for edges
    bool        eyelight    = false;         // camera light mode
    float       exposure    = 0;             // exposure
    float       gamma       = 2.2f;          // gamma
    vec3f       ambient     = {0, 0, 0};     // ambient lighting
    float       near_plane  = 0.01f;         // near plane
    float       far_plane   = 10000.0f;      // far plane

    draw_glstate* state = nullptr;

    bool                                       widgets_open   = false;
    bool                                       navigation_fps = false;
    void*                                      selection      = nullptr;
    std::vector<std::pair<std::string, void*>> update_list;
    float                                      time       = 0;
    std::string                                anim_group = "";
    vec2f                                      time_range = zero2f;
    bool                                       animate    = false;

    ~app_state() {
        if (scn) delete scn;
        if (state) delete state;
    }
};

void draw_glscene(const draw_glstate* state, const scene* scn,
    const camera* cam, const vec2i& viewport_size, const void* highlighted,
    bool eyelight, bool wireframe, bool edges, float exposure, float gamma,
    float near_plane, float far_plane);

// draw with shading
void draw(glwindow* win) {
    auto app              = (app_state*)get_user_pointer(win);
    auto framebuffer_size = get_glframebuffer_size(win);
    app->resolution       = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app->update_list) {
        if (sel.first == "texture") {
            // TODO: update texture
            printf("texture update not supported\n");
        }
        if (sel.first == "subdiv") {
            // TODO: update subdiv
            printf("subdiv update not supported\n");
        }
        if (sel.first == "shape") {
            // TODO: update shape
            printf("shape update not supported\n");
        }
        if (sel.first == "node" || sel.first == "animation" ||
            app->time != last_time) {
            update_transforms(app->scn, app->time, app->anim_group);
            last_time = app->time;
        }
    }
    app->update_list.clear();

    auto cam = app->scn->cameras.at(app->camid);
    clear_glframebuffer(vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    draw_glscene(app->state, app->scn, cam, framebuffer_size, app->selection,
        app->eyelight, app->wireframe, app->edges, app->exposure, app->gamma,
        app->near_plane, app->far_plane);

    begin_glwidgets_frame(win);
    if (begin_glwidgets_window(win, "yview")) {
        if (begin_header_glwidget(win, "scene")) {
            draw_label_glwidgets(win, "scene", "%s", app->filename.c_str());
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "view")) {
            draw_combobox_glwidget(win, "camera", cam, app->scn->cameras, false);
            draw_slider_glwidget(win, "resolution", app->resolution, 256, 4096);
            draw_checkbox_glwidget(win, "eyelight", app->eyelight);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "wireframe", app->wireframe);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "edges", app->edges);
            if (app->time_range != zero2f) {
                draw_slider_glwidget(
                    win, "time", app->time, app->time_range.x, app->time_range.y);
                draw_inputtext_glwidget(win, "anim group", app->anim_group);
                draw_checkbox_glwidget(win, "animate", app->animate);
            }
            draw_slider_glwidget(win, "exposure", app->exposure, -10, 10);
            draw_slider_glwidget(win, "gamma", app->gamma, 0.1f, 4);
            draw_slider_glwidget(win, "near", app->near_plane, 0.01f, 1.0f);
            draw_slider_glwidget(win, "far", app->far_plane, 1000.0f, 10000.0f);
            draw_checkbox_glwidget(win, "fps", app->navigation_fps);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "navigate")) {
            draw_glwidgets_scene_tree(
                win, "", app->scn, app->selection, app->update_list, 200);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "inspect")) {
            draw_glwidgets_scene_inspector(
                win, "", app->scn, app->selection, app->update_list, 200);
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

        void eval_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
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

        bool eval_material(vec2 texcoord, vec4 color, out vec3 ke, 
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
            vec3 txt = 2 * texture(mat_norm_txt,texcoord).xyz - 1;
            txt.y = -txt.y;
            return normalize( tangu * txt.x + tangv * txt.y + norm * txt.z );
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
            bool has_brdf = eval_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

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
                        eval_light(lid, pos, cl, wi);
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
void draw_glshape(const draw_glstate* state, const shape* shp,
    const material* mat, const mat4f& xform, bool highlighted, bool eyelight,
    bool edges) {
    set_gluniform(state->gl_prog, "shape_xform", xform);
    set_gluniform(state->gl_prog, "shape_normal_offset", 0.0f);

    auto mtype = 1;
    if (mat->base_metallic) mtype = 2;
    if (mat->gltf_textures) mtype = (mat->base_metallic) ? 2 : 3;
    set_gluniform(state->gl_prog, "mat_type", mtype);
    set_gluniform(state->gl_prog, "mat_ke", mat->ke);
    set_gluniform(state->gl_prog, "mat_kd", mat->kd);
    set_gluniform(state->gl_prog, "mat_ks", mat->ks);
    set_gluniform(state->gl_prog, "mat_rs", mat->rs);
    set_gluniform(state->gl_prog, "mat_op", mat->op);
    set_gluniform(state->gl_prog, "mat_double_sided", (int)mat->double_sided);
    set_gluniform_texture(state->gl_prog, "mat_ke_txt", "mat_ke_txt_on",
        mat->ke_txt ? state->txt_id.at(mat->ke_txt) : 0, 0);
    set_gluniform_texture(state->gl_prog, "mat_kd_txt", "mat_kd_txt_on",
        mat->kd_txt ? state->txt_id.at(mat->kd_txt) : 0, 1);
    set_gluniform_texture(state->gl_prog, "mat_ks_txt", "mat_ks_txt_on",
        mat->ks_txt ? state->txt_id.at(mat->ks_txt) : 0, 2);
    set_gluniform_texture(state->gl_prog, "mat_rs_txt", "mat_rs_txt_on",
        mat->rs_txt ? state->txt_id.at(mat->rs_txt) : 0, 3);
    set_gluniform_texture(state->gl_prog, "mat_op_txt", "mat_op_txt_on",
        mat->op_txt ? state->txt_id.at(mat->op_txt) : 0, 4);
    set_gluniform_texture(state->gl_prog, "mat_norm_txt", "mat_norm_txt_on",
        mat->norm_txt ? state->txt_id.at(mat->norm_txt) : 0, 5);

    auto& vbos = state->shp_vbos.at(shp);
    set_gluniform(state->gl_prog, "elem_faceted", (int)shp->norm.empty());
    set_glvertexattrib(state->gl_prog, "vert_pos", vbos.gl_pos, zero3f);
    set_glvertexattrib(state->gl_prog, "vert_norm", vbos.gl_norm, zero3f);
    set_glvertexattrib(
        state->gl_prog, "vert_texcoord", vbos.gl_texcoord, zero2f);
    set_glvertexattrib(
        state->gl_prog, "vert_color", vbos.gl_color, vec4f{1, 1, 1, 1});
    set_glvertexattrib(
        state->gl_prog, "vert_tangsp", vbos.gl_tangsp, vec4f{0, 0, 1, 1});

    if (!shp->points.empty()) {
        set_gluniform(state->gl_prog, "elem_type", 1);
        draw_glpoints(vbos.gl_points, shp->points.size());
    }
    if (!shp->lines.empty()) {
        set_gluniform(state->gl_prog, "elem_type", 2);
        draw_glpoints(vbos.gl_lines, shp->lines.size());
    }
    if (!shp->triangles.empty()) {
        set_gluniform(state->gl_prog, "elem_type", 3);
        draw_gltriangles(vbos.gl_triangles, shp->triangles.size());
    }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_glerror();
        set_gluniform(state->gl_prog, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state->gl_prog, "ke"), 0, 0, 0);
        set_gluniform(state->gl_prog, "op"), mat->op);
        set_gluniform(state->gl_prog, "shp_normal_offset"), 0.01f);
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

// Display a scene
void draw_glscene(const draw_glstate* state, const scene* scn,
    const camera* cam, const vec2i& viewport_size, const void* highlighted,
    bool eyelight, bool wireframe, bool edges, float exposure, float gamma,
    float near_plane, float far_plane) {
    set_glviewport(viewport_size);

    auto camera_view = frame_to_mat(inverse(cam->frame));
    // auto camera_proj =
    //         perspective_mat(eval_camera_fovy(cam),
    //             (float)viewport_size.x / (float)viewport_size.y, near_plane);
    auto camera_proj = perspective_mat(eval_camera_fovy(cam),
        (float)viewport_size.x / (float)viewport_size.y, near_plane, far_plane);

    bind_glprogram(state->gl_prog);
    set_gluniform(state->gl_prog, "cam_pos", cam->frame.o);
    set_gluniform(state->gl_prog, "cam_xform_inv", camera_view);
    set_gluniform(state->gl_prog, "cam_proj", camera_proj);
    set_gluniform(state->gl_prog, "eyelight", (int)eyelight);
    set_gluniform(state->gl_prog, "exposure", exposure);
    set_gluniform(state->gl_prog, "gamma", gamma);

    if (!eyelight) {
        auto lights_pos  = std::vector<vec3f>();
        auto lights_ke   = std::vector<vec3f>();
        auto lights_type = std::vector<int>();
        for (auto lgt : scn->instances) {
            if (lgt->mat->ke == zero3f) continue;
            if (lights_pos.size() >= 16) break;
            auto shp  = lgt->shp;
            auto bbox = compute_bbox(shp);
            auto pos  = (bbox.max + bbox.min) / 2;
            auto area = 0.0f;
            if (!shp->triangles.empty()) {
                for (auto t : shp->triangles)
                    area += triangle_area(
                        shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
            } else if (!shp->lines.empty()) {
                for (auto l : shp->lines)
                    area += line_length(shp->pos[l.x], shp->pos[l.y]);
            } else {
                area += shp->pos.size();
            }
            auto ke = lgt->mat->ke * area;
            lights_pos.push_back(transform_point(lgt->frame, pos));
            lights_ke.push_back(ke);
            lights_type.push_back(0);
        }
        if (lights_pos.empty()) eyelight = false;
        set_gluniform(state->gl_prog, "lamb", zero3f);
        set_gluniform(state->gl_prog, "lnum", (int)lights_pos.size());
        for (auto i = 0; i < lights_pos.size(); i++) {
            auto is = std::to_string(i);
            set_gluniform(
                state->gl_prog, ("lpos[" + is + "]").c_str(), lights_pos[i]);
            set_gluniform(
                state->gl_prog, ("lke[" + is + "]").c_str(), lights_ke[i]);
            set_gluniform(state->gl_prog, ("ltype[" + is + "]").c_str(),
                (int)lights_type[i]);
        }
    }

    if (wireframe) set_glwireframe(true);
    for (auto ist : scn->instances) {
        draw_glshape(state, ist->shp, ist->mat, frame_to_mat(ist->frame),
            ist == highlighted || ist->shp == highlighted ||
                ist->mat == highlighted,
            eyelight, edges);
    }

    bind_glprogram(0);
    if (wireframe) set_glwireframe(false);
}

draw_glstate* init_draw_state(glwindow* win) {
    auto app   = (app_state*)get_user_pointer(win);
    auto state = new draw_glstate();
    // load textures and vbos
    state->gl_prog         = make_glprogram(vertex, fragment);
    state->txt_id[nullptr] = 0;
    for (auto txt : app->scn->textures) {
        if (!txt->imgf.empty()) {
            state->txt_id[txt] = make_gltexture(txt->imgf, true, true, true);
        } else if (!txt->imgb.empty()) {
            state->txt_id[txt] = make_gltexture(txt->imgb, txt->srgb, true, true);
        } else {
            printf("bad texture");
        }
    }
    for (auto& shp : app->scn->shapes) {
        auto vbos = draw_glshape_vbos();
        if (!shp->pos.empty())
            vbos.gl_pos = make_glarraybuffer(shp->pos, false);
        if (!shp->norm.empty())
            vbos.gl_norm = make_glarraybuffer(shp->norm, false);
        if (!shp->texcoord.empty())
            vbos.gl_texcoord = make_glarraybuffer(shp->texcoord, false);
        if (!shp->color.empty())
            vbos.gl_color = make_glarraybuffer(shp->color, false);
        if (!shp->tangsp.empty())
            vbos.gl_tangsp = make_glarraybuffer(shp->tangsp, false);
        if (!shp->points.empty())
            vbos.gl_points = make_glelementbuffer(shp->points, false);
        if (!shp->lines.empty())
            vbos.gl_lines = make_glelementbuffer(shp->lines, false);
        if (!shp->triangles.empty())
            vbos.gl_triangles = make_glelementbuffer(shp->triangles, false);
        state->shp_vbos[shp] = vbos;
    }
    return state;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto cam   = app->scn->cameras.at(app->camid);
    auto wsize = clamp(eval_image_size(cam, app->resolution), 256, 1440);
    auto win   = make_glwindow(wsize, "yview", app, draw);

    // init widget
    init_glwidgets(win);

    // load textures and vbos
    update_transforms(app->scn, app->time);

    // init gl data
    app->state = init_draw_state(win);

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
            auto cam = app->scn->cameras.at(app->camid);
            camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back({"camera", cam});
        }

        // animation
        if (app->animate) {
            app->time += 1 / 60.0f;
            if (app->time < app->time_range.x || app->time > app->time_range.y)
                app->time = app->time_range.x;
            update_transforms(app->scn, app->time);
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
std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
load_ini(const std::string& filename) {
    auto f = fopen(filename.c_str(), "rt");
    if (!f) throw std::runtime_error("cannot open " + filename);
    auto ret       = std::unordered_map<std::string,
        std::unordered_map<std::string, std::string>>();
    auto cur_group = std::string();
    ret[""]        = {};

    char buf[4096];
    while (fgets(buf, 4096, f)) {
        auto line = std::string(buf);
        if (line.empty()) continue;
        if (line.front() == ';') continue;
        if (line.front() == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') throw std::runtime_error("bad INI format");
            cur_group      = line.substr(1, line.length() - 2);
            ret[cur_group] = {};
        } else if (line.find('=') != line.npos) {
            auto var            = line.substr(0, line.find('='));
            auto val            = line.substr(line.find('=') + 1);
            ret[cur_group][var] = val;
        } else {
            throw std::runtime_error("bad INI format");
        }
    }

    fclose(f);

    return ret;
}

int main(int argc, char* argv[]) {
    // initialize app
    auto app = new app_state();

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "views scenes inteactively", "yview");
    app->camid      = parse_arg(parser, "--camera", 0, "Camera index.");
    app->resolution = parse_arg(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->eyelight = parse_arg(
        parser, "--eyelight,-c", false, "Eyelight rendering.");
    auto double_sided = parse_arg(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto quiet = parse_arg(
        parser, "--quiet,-q", false, "Print only errors messages");
    auto highlight_filename = parse_arg(
        parser, "--highlights", "", "Highlight filename");
    app->imfilename = parse_arg(
        parser, "--output-image,-o", "out.png", "Image filename");
    app->filename = parse_arg(
        parser, "scene", "scene.json", "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    if (!quiet) printf("loading scene %s\n", app->filename.c_str());
    try {
        app->scn = load_scene(app->filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\n", app->filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    tesselate_subdivs(app->scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (double_sided) {
        for (auto mat : app->scn->materials) mat->double_sided = true;
    }
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            make_bbox_camera("<view>", compute_bbox(app->scn)));
    add_missing_names(app->scn);
    for (auto& err : validate(app->scn)) printf("warning: %s\n", err.c_str());

    // animation
    auto time_range = compute_animation_range(app->scn);
    app->time       = time_range.x;

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
