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

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#define GLFW_INCLUDE_GLCOREARB
#else
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#endif
#include <GLFW/glfw3.h>

#include "ysceneui.h"

#include "imgui/imgui.h"
#include "imgui/imgui_ext.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

// Application state
struct app_state {
    // scene
    ygl::scene* scn = nullptr;

    // parameters
    std::string filename = "scene.json";            // scene name
    std::string imfilename = "out.png";             // output image
    std::string outfilename = "scene.json";         // save scene name
    int camid = 0;                                  // camera id
    int resolution = 512;                           // image resolution
    bool wireframe = false;                         // wireframe drawing
    bool edges = false;                             // draw edges
    float edge_offset = 0.01f;                      // offset for edges
    bool eyelight = false;                          // camera light mode
    float exposure = 0;                             // exposure
    float gamma = 2.2f;                             // gamma
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};  // background
    ygl::vec3f ambient = {0, 0, 0};                 // ambient lighting

    unsigned int gl_prog = 0;

    bool widgets_open = false;
    bool navigation_fps = false;
    void* selection = nullptr;
    std::vector<std::pair<std::string, void*>> update_list;
    float time = 0;
    std::string anim_group = "";
    ygl::vec2f time_range = ygl::zero2f;
    bool animate = false;

    ~app_state() {
        if (scn) delete scn;
    }
};

void draw_glscene(const ygl::scene* scn, const ygl::camera* cam,
    unsigned int prog, const ygl::vec2i& viewport_size, const void* highlighted,
    bool eyelight, bool wireframe, bool edges, float exposure, float gamma);

// draw with shading
void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    auto window_size = ygl::zero2i;
    auto framebuffer_size = ygl::zero2i;
    glfwGetWindowSize(win, &window_size.x, &window_size.y);
    glfwGetFramebufferSize(win, &framebuffer_size.x, &framebuffer_size.y);
    app->resolution = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app->update_list) {
        if (sel.first == "texture") {
            // TODO: update texture
            printf("texture update not supported\n");
        }
        if (sel.first == "subdiv") {
            // TODO: update subdiv
            printf("texture update not supported\n");
        }
        if (sel.first == "shape") {
            // TODO: update shape
            printf("texture update not supported\n");
        }
        if (sel.first == "node" || sel.first == "animation" ||
            app->time != last_time) {
            ygl::update_transforms(app->scn, app->time, app->anim_group);
            last_time = app->time;
        }
        if (sel.first == "shape" || sel.first == "material" ||
            sel.first == "node") {
            ygl::init_lights(app->scn);
            if (app->scn->lights.empty()) app->eyelight = true;
        }
    }
    app->update_list.clear();

    auto cam = app->scn->cameras.at(app->camid);
    glClearColor(app->background.x, app->background.y, app->background.z,
        app->background.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    draw_glscene(app->scn, cam, app->gl_prog, framebuffer_size, app->selection,
        app->eyelight, app->wireframe, app->edges, app->exposure, app->gamma);

    static auto first_time = true;
    // auto app = (app_state*)glfwGetWindowUserPointer(win);
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }

    if (ImGui::Begin("yview")) {
        ImGui::LabelText("scene", "%s", app->filename.c_str());
        if (app->time_range != ygl::zero2f) {
            ImGui::SliderFloat(
                "time", &app->time, app->time_range.x, app->time_range.y);
            ImGui::InputText("anim group", &app->anim_group);
            ImGui::Checkbox("animate", &app->animate);
        }
        if (ImGui::TreeNode("render settings")) {
            ImGui::Combo("camera", &cam, app->scn->cameras, false);
            ImGui::SliderInt("resolution", &app->resolution, 256, 4096);
            ImGui::Checkbox("eyelight", &app->eyelight);
            ImGui::SameLine();
            ImGui::Checkbox("wireframe", &app->wireframe);
            ImGui::SameLine();
            ImGui::Checkbox("edges", &app->edges);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("view settings")) {
            ImGui::SliderFloat("exposure", &app->exposure, -10, 10);
            ImGui::SliderFloat("gamma", &app->gamma, 0.1f, 4);
            ImGui::ColorEdit4("background", &app->background.x);
            ImGui::Checkbox("fps", &app->navigation_fps);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene tree")) {
            draw_glwidgets_scene_tree(
                "", app->scn, app->selection, app->update_list, 200);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene object")) {
            draw_glwidgets_scene_inspector(
                "", app->scn, app->selection, app->update_list, 200);
            ImGui::TreePop();
        }
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(win);
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

        float pi = 3.14159265;

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
                return ((1+dot(wo,wi))/2) * kd/pi;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pi;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pi)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pi;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pi)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
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

            ke_txt.xyz = pow(ke_txt.xyz, vec3(2.2));
            kd_txt.xyz = pow(kd_txt.xyz, vec3(2.2));
            ks_txt.xyz = pow(ks_txt.xyz, vec3(2.2));

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
                    c += pi * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
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
void draw_glshape(const ygl::shape* shp, const ygl::material* mat,
    const ygl::mat4f& xform, bool highlighted, unsigned int prog, bool eyelight,
    bool edges) {
    glUniformMatrix4fv(
        glGetUniformLocation(prog, "shape_xform"), 1, false, &xform.x.x);
    glUniform1f(glGetUniformLocation(prog, "shape_normal_offset"), 0.0f);

    auto uniform_texture = [](auto& prog, const char* name, const char* name_on,
                               const ygl::texture* txt, int unit) {
        if (txt) {
            glActiveTexture(GL_TEXTURE0 + unit);
            glBindTexture(GL_TEXTURE_2D, txt->gl_txt);
            glUniform1i(glGetUniformLocation(prog, name), unit);
            glUniform1i(glGetUniformLocation(prog, name_on), 1);
        } else {
            glUniform1i(glGetUniformLocation(prog, name_on), 0);
        }
    };

    auto mtype = 1;
    if (mat->base_metallic) mtype = 2;
    if (mat->gltf_textures) mtype = (mat->base_metallic) ? 2 : 3;
    glUniform1i(glGetUniformLocation(prog, "mat_type"), mtype);
    glUniform3f(
        glGetUniformLocation(prog, "mat_ke"), mat->ke.x, mat->ke.y, mat->ke.z);
    glUniform3f(
        glGetUniformLocation(prog, "mat_kd"), mat->kd.x, mat->kd.y, mat->kd.z);
    glUniform3f(
        glGetUniformLocation(prog, "mat_ks"), mat->ks.x, mat->ks.y, mat->ks.z);
    glUniform1f(glGetUniformLocation(prog, "mat_rs"), mat->rs);
    glUniform1f(glGetUniformLocation(prog, "mat_op"), mat->op);
    glUniform1i(
        glGetUniformLocation(prog, "mat_double_sided"), (int)mat->double_sided);
    uniform_texture(prog, "mat_ke_txt", "mat_ke_txt_on", mat->ke_txt, 0);
    uniform_texture(prog, "mat_kd_txt", "mat_kd_txt_on", mat->kd_txt, 1);
    uniform_texture(prog, "mat_ks_txt", "mat_ks_txt_on", mat->ks_txt, 2);
    uniform_texture(prog, "mat_rs_txt", "mat_rs_txt_on", mat->rs_txt, 3);
    uniform_texture(prog, "mat_op_txt", "mat_op_txt_on", mat->op_txt, 4);
    uniform_texture(prog, "mat_norm_txt", "mat_norm_txt_on", mat->norm_txt, 5);

    glUniform1i(
        glGetUniformLocation(prog, "elem_faceted"), (int)shp->norm.empty());
    if (shp->gl_pos) {
        glEnableVertexAttribArray(glGetAttribLocation(prog, "vert_pos"));
        glBindBuffer(GL_ARRAY_BUFFER, shp->gl_pos);
        glVertexAttribPointer(
            glGetAttribLocation(prog, "vert_pos"), 3, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib3f(glGetAttribLocation(prog, "vert_pos"), 0, 0, 0);
    }
    if (shp->gl_norm) {
        glEnableVertexAttribArray(glGetAttribLocation(prog, "vert_norm"));
        glBindBuffer(GL_ARRAY_BUFFER, shp->gl_norm);
        glVertexAttribPointer(
            glGetAttribLocation(prog, "vert_norm"), 3, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib3f(glGetAttribLocation(prog, "vert_norm"), 0, 0, 0);
    }
    if (shp->gl_texcoord) {
        glEnableVertexAttribArray(glGetAttribLocation(prog, "vert_texcoord"));
        glBindBuffer(GL_ARRAY_BUFFER, shp->gl_texcoord);
        glVertexAttribPointer(glGetAttribLocation(prog, "vert_texcoord"), 2,
            GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib2f(glGetAttribLocation(prog, "vert_texcoord"), 0, 0);
    }
    if (shp->gl_color) {
        glEnableVertexAttribArray(glGetAttribLocation(prog, "vert_color"));
        glBindBuffer(GL_ARRAY_BUFFER, shp->gl_color);
        glVertexAttribPointer(
            glGetAttribLocation(prog, "vert_color"), 4, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib4f(glGetAttribLocation(prog, "vert_color"), 1, 1, 1, 1);
    }
    if (shp->gl_tangsp) {
        glEnableVertexAttribArray(glGetAttribLocation(prog, "vert_tangsp"));
        glBindBuffer(GL_ARRAY_BUFFER, shp->gl_tangsp);
        glVertexAttribPointer(
            glGetAttribLocation(prog, "vert_tangsp"), 4, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib4f(glGetAttribLocation(prog, "vert_tangsp"), 0, 0, 1, 1);
    }

    if (!shp->points.empty()) {
        glUniform1i(glGetUniformLocation(prog, "elem_type"), 1);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shp->gl_points);
        glDrawElements(GL_POINTS, 1 * shp->points.size(), GL_UNSIGNED_INT, 0);
    }
    if (!shp->lines.empty()) {
        glUniform1i(glGetUniformLocation(prog, "elem_type"), 2);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shp->gl_lines);
        glDrawElements(GL_LINES, 2 * shp->lines.size(), GL_UNSIGNED_INT, 0);
    }
    if (!shp->triangles.empty()) {
        glUniform1i(glGetUniformLocation(prog, "elem_type"), 3);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shp->gl_triangles);
        glDrawElements(
            GL_TRIANGLES, 3 * shp->triangles.size(), GL_UNSIGNED_INT, 0);
    }

#if 0
    if ((shp->gl_edges && edges && !wireframe) || highlighted) {
        ygl::enable_glculling(false);
        check_glerror();
        glUniform1i(glGetUniformLocation(prog, "mtype"), 0);
        glUniform3f(glGetUniformLocation(prog, "ke"), 0, 0, 0);
        glUniform1f(glGetUniformLocation(prog, "op"), mat->op);
        glUniform1f(glGetUniformLocation(prog, "shp_normal_offset"), 0.01f);
        check_glerror();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shp->gl_edges);
        glDrawElements(GL_LINES, shp->triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_glerror();
    }
#endif
    if (edges) printf("edges are momentarily disabled\n");

    for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
}

// Display a scene
void draw_glscene(const ygl::scene* scn, const ygl::camera* cam,
    unsigned int prog, const ygl::vec2i& viewport_size, const void* highlighted,
    bool eyelight, bool wireframe, bool edges, float exposure, float gamma) {
    glViewport(0, 0, viewport_size.x, viewport_size.y);

    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj =
        (cam->far >= ygl::flt_max) ?
            ygl::perspective_mat(eval_camera_fovy(cam),
                (float)viewport_size.x / (float)viewport_size.y, cam->near) :
            ygl::perspective_mat(eval_camera_fovy(cam),
                (float)viewport_size.x / (float)viewport_size.y, cam->near,
                cam->far);

    glUseProgram(prog);
    glUniform3f(glGetUniformLocation(prog, "cam_pos"), cam->frame.o.x,
        cam->frame.o.y, cam->frame.o.z);
    glUniformMatrix4fv(glGetUniformLocation(prog, "cam_xform_inv"), 1, false,
        &camera_view.x.x);
    glUniformMatrix4fv(
        glGetUniformLocation(prog, "cam_proj"), 1, false, &camera_proj.x.x);
    glUniform1i(glGetUniformLocation(prog, "eyelight"), (int)eyelight);
    glUniform1f(glGetUniformLocation(prog, "exposure"), exposure);
    glUniform1f(glGetUniformLocation(prog, "gamma"), gamma);

    if (!eyelight) {
        auto lights_pos = std::vector<ygl::vec3f>();
        auto lights_ke = std::vector<ygl::vec3f>();
        auto lights_type = std::vector<int>();
        for (auto lgt : scn->lights) {
            if (lights_pos.size() >= 16) break;
            auto shp = lgt->shp;
            auto bbox = compute_bbox(shp);
            auto pos = (bbox.max + bbox.min) / 2;
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
        glUniform3f(glGetUniformLocation(prog, "lamb"), 0, 0, 0);
        glUniform1i(glGetUniformLocation(prog, "lnum"), (int)lights_pos.size());
        for (auto i = 0; i < lights_pos.size(); i++) {
            auto is = std::to_string(i);
            glUniform3f(
                glGetUniformLocation(prog, ("lpos[" + is + "]").c_str()),
                lights_pos[i].x, lights_pos[i].y, lights_pos[i].z);
            glUniform3f(glGetUniformLocation(prog, ("lke[" + is + "]").c_str()),
                lights_ke[i].x, lights_ke[i].y, lights_ke[i].z);
            glUniform1i(
                glGetUniformLocation(prog, ("ltype[" + is + "]").c_str()),
                (int)lights_type[i]);
        }
    }

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (auto ist : scn->instances) {
        draw_glshape(ist->shp, ist->mat, frame_to_mat(ist->frame),
            ist == highlighted || ist->shp == highlighted ||
                ist->mat == highlighted,
            prog, eyelight, edges);
    }

    glUseProgram(0);
    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

inline unsigned int make_glprogram(const char* vertex, const char* fragment,
    unsigned int& vid, unsigned int& fid, unsigned int& vao) {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    int errflags;
    char errbuf[10000];

    // create vertex
    vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vid, 1, &vertex, NULL);
    glCompileShader(vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }

    // create fragment
    fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fid, 1, &fragment, NULL);
    glCompileShader(fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }

    // create program
    auto pid = glCreateProgram();
    glAttachShader(pid, vid);
    glAttachShader(pid, fid);
    glLinkProgram(pid);
    glValidateProgram(pid);
    glGetProgramiv(pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }

    return pid;
}

inline unsigned int make_glprogram(const char* vertex, const char* fragment) {
    unsigned int vid = 0, fid = 0, vao = 0;
    return make_glprogram(vertex, fragment, vid, fid, vao);
}

void init_drawscene(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    // load textures and vbos
    app->gl_prog = make_glprogram(vertex, fragment);
    for (auto& txt : app->scn->textures) {
        glGenTextures(1, &txt->gl_txt);
        glBindTexture(GL_TEXTURE_2D, txt->gl_txt);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txt->img.width, txt->img.height,
            0, GL_RGBA, GL_FLOAT, txt->img.pxl.data());
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    auto make_glbuffer = [](auto& data, bool elems) {
        auto bid = (unsigned int)0;
        auto target = (elems) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
        glGenBuffers(1, &bid);
        glBindBuffer(target, bid);
        glBufferData(target, sizeof(data.at(0)) * data.size(), data.data(),
            GL_DYNAMIC_DRAW);
        glBindBuffer(target, 0);
        return bid;
    };
    for (auto& shp : app->scn->shapes) {
        if (!shp->pos.empty()) shp->gl_pos = make_glbuffer(shp->pos, false);
        if (!shp->norm.empty()) shp->gl_norm = make_glbuffer(shp->norm, false);
        if (!shp->texcoord.empty())
            shp->gl_texcoord = make_glbuffer(shp->texcoord, false);
        if (!shp->color.empty())
            shp->gl_color = make_glbuffer(shp->color, false);
        if (!shp->tangsp.empty())
            shp->gl_tangsp = make_glbuffer(shp->tangsp, false);
        if (!shp->points.empty())
            shp->gl_points = make_glbuffer(shp->points, true);
        if (!shp->lines.empty())
            shp->gl_lines = make_glbuffer(shp->lines, true);
        if (!shp->triangles.empty())
            shp->gl_triangles = make_glbuffer(shp->triangles, true);
    }
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto cam = app->scn->cameras.at(app->camid);
    auto ww = ygl::clamp(ygl::image_width(cam, app->resolution), 256, 1440);
    auto wh = ygl::clamp(app->resolution, 256, 1440);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    auto win = glfwCreateWindow(ww, wh, "yview", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    // init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif

    glfwSetWindowRefreshCallback(win, draw);
    glfwSetWindowUserPointer(win, app);

    // init widget
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename = nullptr;
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init();
    ImGui::StyleColorsDark();

    // load textures and vbos
    ygl::update_transforms(app->scn, app->time);
    ygl::init_lights(app->scn);
    if (app->scn->lights.empty()) app->eyelight = true;

    // init gl data
    init_drawscene(win);

    // loop
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        double mouse_posx, mouse_posy;
        glfwGetCursorPos(win, &mouse_posx, &mouse_posy);
        mouse_pos = {(float)mouse_posx, (float)mouse_posy};
        auto mouse_left =
            glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        auto mouse_right =
            glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
        auto alt_down = glfwGetKey(win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
                        glfwGetKey(win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
        auto shift_down = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse and keyboard for navigation
        if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->camid);
            ygl::camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back({"camera", cam});
        }

        // animation
        if (app->animate) {
            app->time += 1 / 60.0f;
            if (app->time < app->time_range.x || app->time > app->time_range.y)
                app->time = app->time_range.x;
            ygl::update_transforms(app->scn, app->time);
        }

        // draw
        draw(win);

        // event hadling
        if ((mouse_left || mouse_right) || widgets_active) {
            glfwPollEvents();
        } else {
            glfwWaitEvents();
        }
    }
}

// Load INI file. The implementation does not handle escaping.
std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
load_ini(const std::string& filename) {
    auto f = fopen(filename.c_str(), "rt");
    if (!f) throw std::runtime_error("cannot open " + filename);
    auto ret = std::unordered_map<std::string,
        std::unordered_map<std::string, std::string>>();
    auto cur_group = std::string();
    ret[""] = {};

    char buf[4096];
    while (fgets(buf, 4096, f)) {
        auto line = std::string(buf);
        if (line.empty()) continue;
        if (line.front() == ';') continue;
        if (line.front() == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') throw std::runtime_error("bad INI format");
            cur_group = line.substr(1, line.length() - 2);
            ret[cur_group] = {};
        } else if (line.find('=') != line.npos) {
            auto var = line.substr(0, line.find('='));
            auto val = line.substr(line.find('=') + 1);
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
    auto parser = ygl::make_cmdline_parser(
        argc, argv, "views scenes inteactively", "yview");
    app->camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    app->resolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->eyelight =
        ygl::parse_flag(parser, "--eyelight,-c", false, "Eyelight rendering.");
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto quiet = ygl::parse_flag(
        parser, "--quiet,-q", false, "Print only errors messages");
    auto highlight_filename =
        ygl::parse_string(parser, "--highlights", "", "Highlight filename");
    app->imfilename = ygl::parse_string(
        parser, "--output-image,-o", "out.png", "Image filename");
    app->filename = ygl::parse_string(
        parser, "scene", "scene.json", "Scene filename", true);
    ygl::check_cmdline(parser);

    // scene loading
    if (!quiet) printf("loading scene %s\n", app->filename.c_str());
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\n", app->filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(app->scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (double_sided) {
        for (auto mat : app->scn->materials) mat->double_sided = true;
    }
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", ygl::compute_bbox(app->scn)));
    ygl::add_missing_names(app->scn);
    for (auto& err : ygl::validate(app->scn))
        printf("warning: %s\n", err.c_str());

    // animation
    auto time_range = ygl::compute_animation_range(app->scn);
    app->time = time_range.x;

    // lights
    ygl::init_lights(app->scn);

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
