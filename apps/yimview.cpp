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
#include "CLI11.hpp"
#include "yglui.h"
using namespace std::literals;

#include <memory>
#include <unordered_map>

// Generic image that contains either an HDR or an LDR image, giving access
// to both. This is helpful when writing viewers or generic image
// manipulation code
struct gimage {
    std::string name = ""s;
    std::string filename = ""s;  // path
    ygl::image4f img = {};       // image

    // min/max values
    ygl::bbox4f pxl_bounds = ygl::invalid_bbox4f;
    ygl::bbox1f lum_bounds = ygl::invalid_bbox1f;

    // image adjustment
    bool updated = true;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    bool is_hdr = false;

    // display image
    ygl::image4f display;

    // opengl texture
    uint gl_txt = 0;
};

struct app_state {
    std::vector<std::shared_ptr<gimage>> imgs;
    std::shared_ptr<gimage> img = nullptr;

    // viewing properties
    uint gl_prog = 0, gl_vbo = 0, gl_ebo = 0;
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
};

// compute min/max
void update_minmax(std::shared_ptr<gimage> img) {
    img->pxl_bounds = ygl::invalid_bbox4f;
    img->lum_bounds = ygl::invalid_bbox1f;
    for (auto& p : img->img) {
        img->pxl_bounds += p;
        img->lum_bounds += ygl::luminance(xyz(p));
    }
}

// Loads a generic image
std::shared_ptr<gimage> load_gimage(
    const std::string& filename, float exposure, float gamma, bool filmic) {
    auto img = std::make_shared<gimage>();
    img->filename = filename;
    img->name = ygl::get_filename(filename);
    img->exposure = exposure;
    img->gamma = gamma;
    img->filmic = filmic;
    img->is_hdr = ygl::is_hdr_filename(filename);
    img->img = ygl::load_image(filename);
    update_minmax(img);
    return img;
}

std::shared_ptr<gimage> diff_gimage(
    std::shared_ptr<gimage> a, std::shared_ptr<gimage> b, bool color) {
    if (a->img.width() != b->img.width() || a->img.height() != b->img.height())
        return nullptr;
    auto d = std::make_shared<gimage>();
    d->name = "diff " + a->name + " " + b->name;
    d->filename = "";
    d->img = a->img;
    if (color) {
        for (auto i = 0; i < a->img.size(); i++) {
            d->img[i] = {std::abs(a->img[i].x - b->img[i].x),
                std::abs(a->img[i].y - b->img[i].y),
                std::abs(a->img[i].z - b->img[i].z),
                std::max(a->img[i].w, b->img[i].w)};
        }
    } else {
        for (auto i = 0; i < a->img.size(); i++) {
            auto la = (a->img[i].x + a->img[i].y + a->img[i].z) / 3;
            auto lb = (b->img[i].x + b->img[i].y + b->img[i].z) / 3;
            auto ld = fabsf(la - lb);
            d->img[i] = {ld, ld, ld, std::max(a->img[i].w, b->img[i].w)};
        }
    }
    update_minmax(d);
    return d;
}

void update_display_image(const std::shared_ptr<gimage>& img) {
    img->display = img->img;
    if (img->is_hdr) {
        img->display = ygl::tonemap_image(
            img->display, img->exposure, img->gamma, img->filmic);
    }
}

void draw_widgets(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    if (begin_widgets_frame(win, "yimview")) {
        ImGui::Combo("image", &app->img, app->imgs, false);
        ImGui::LabelText("filename", "%s", app->img->filename.c_str());
        ImGui::LabelText(
            "size", "%d x %d ", app->img->img.width(), app->img->img.height());
        auto edited = 0;
        edited += ImGui::SliderFloat("exposure", &app->img->exposure, -5, 5);
        edited += ImGui::SliderFloat("gamma", &app->img->gamma, 1, 3);
        edited += ImGui::Checkbox("filmic", &app->img->filmic);
        if (edited) app->img->updated = true;
        auto zoom = app->imframe.x.x;
        if (ImGui::SliderFloat("zoom", &zoom, 0.1, 10))
            app->imframe.x.x = app->imframe.y.y = zoom;
        ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
        ImGui::ColorEdit4("background", &app->background.x);
        ImGui::Separator();
        auto mouse_x = 0.0, mouse_y = 0.0;
        glfwGetCursorPos(win, &mouse_x, &mouse_y);
        auto ij =
            ygl::get_image_coords(ygl::vec2f{(float)mouse_x, (float)mouse_y},
                app->imframe, {app->img->img.width(), app->img->img.height()});
        ImGui::DragInt2("mouse", &ij.x);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < app->img->img.width() && ij.y >= 0 &&
            ij.y < app->img->img.height()) {
            pixel = app->img->img.at(ij);
        }
        ImGui::ColorEdit4("pixel", &pixel.x);
        if (!app->img->img.empty()) {
            ImGui::DragFloat4("pxl min", &app->img->pxl_bounds.min.x);
            ImGui::DragFloat4("pxl max", &app->img->pxl_bounds.max.y);
            ImGui::DragFloat("lum min", &app->img->lum_bounds.min);
            ImGui::DragFloat("lum max", &app->img->lum_bounds.max);
        }
    }
    end_widgets_frame();
}

void draw_image(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    auto window_size = ygl::zero2i, framebuffer_size = ygl::zero2i;
    glfwGetWindowSize(win, &window_size.x, &window_size.y);
    glfwGetFramebufferSize(win, &framebuffer_size.x, &framebuffer_size.y);

    ygl::center_image(app->imframe,
        {app->img->img.width(), app->img->img.height()}, window_size,
        app->zoom_to_fit);

    assert(glGetError() == GL_NO_ERROR);
    glViewport(0, 0, framebuffer_size.x, framebuffer_size.y);
    glClearColor(app->background.x, app->background.y, app->background.z,
        app->background.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);
    glUseProgram(app->gl_prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, app->img->gl_txt);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glUniform2f(glGetUniformLocation(app->gl_prog, "win_size"), window_size.x,
        window_size.y);
    glUniform2f(glGetUniformLocation(app->gl_prog, "txt_size"),
        app->img->img.width(), app->img->img.height());
    glUniformMatrix3x2fv(glGetUniformLocation(app->gl_prog, "frame"), 1, false,
        &app->imframe.x.x);
    glUniform1i(glGetUniformLocation(app->gl_prog, "img"), 0);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glEnableVertexAttribArray(
        glGetAttribLocation(app->gl_prog, "vert_texcoord"));
    glBindBuffer(GL_ARRAY_BUFFER, app->gl_vbo);
    glVertexAttribPointer(glGetAttribLocation(app->gl_prog, "vert_texcoord"), 2,
        GL_FLOAT, false, 0, 0);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, app->gl_ebo);
    glDrawElements(GL_TRIANGLES, 3 * 4, GL_UNSIGNED_INT, 0);
    assert(glGetError() == GL_NO_ERROR);

    glUseProgram(0);
    glDisable(GL_BLEND);
}

void draw(GLFWwindow* win) {
    draw_image(win);
    draw_widgets(win);
    glfwSwapBuffers(win);
}

static const char* vertex =
    R"(
    #version 330

    layout(location = 0) in vec2 vert_texcoord;

    uniform mat3x2 frame;
    uniform vec2 txt_size;
    uniform vec2 win_size;
    uniform sampler2D img;

    out vec2 texcoord;

    void main() {
        texcoord = vert_texcoord.xy;
        vec2 pos = frame * vec3(txt_size.x * (vert_texcoord.x - 0.5), 
                                txt_size.y * (vert_texcoord.y - 0.5), 1);
        vec2 upos = 2 * pos / win_size - vec2(1,1);
        upos.y = - upos.y;
        gl_Position = vec4(upos.x, upos.y, 0, 1);
    }

    )";

static const char* fragment =
    R"(
    #version 330

    in vec2 texcoord;
    uniform sampler2D img;
    out vec4 color;

    void main() {
        color = texture(img,texcoord);
    }
    )";

void init_drawimage(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    app->gl_prog = make_glprogram(vertex, fragment);
    auto uv = std::vector<ygl::vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    auto triangles = std::vector<ygl::vec3i>{{0, 1, 2}, {0, 2, 3}};
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &app->gl_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, app->gl_vbo);
    glBufferData(GL_ARRAY_BUFFER, 4 * 2 * 4, uv.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glGenBuffers(1, &app->gl_ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, app->gl_ebo);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER, 2 * 3 * 4, triangles.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(glGetError() == GL_NO_ERROR);
    for (auto& img : app->imgs) {
        assert(glGetError() == GL_NO_ERROR);
        glGenTextures(1, &img->gl_txt);
        glBindTexture(GL_TEXTURE_2D, img->gl_txt);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img->img.width(),
            img->img.height(), 0, GL_RGBA, GL_FLOAT, img->img.data());
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
        assert(glGetError() == GL_NO_ERROR);
    }
}

void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win_width = ygl::clamp(app->imgs[0]->img.width(), 512, 1024);
    auto win_height = ygl::clamp(app->imgs[0]->img.height(), 512, 1024);
    auto win = make_window(win_width, win_height, "yimview", app.get(), draw);

    // init widgets
    init_widgets(win);

    // load textures
    init_drawimage(win);

    // window values
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    auto mouse_button = 0;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        glfwGetCursorPosExt(win, &mouse_pos.x, &mouse_pos.y);
        mouse_button = glfwGetMouseButtonIndexExt(win);
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse
        if (mouse_button && !widgets_active) {
            if (mouse_button == 1) app->imframe.o += mouse_pos - last_pos;
            if (mouse_button == 2) {
                auto zoom = powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
                app->imframe =
                    app->imframe * ygl::frame2f{{zoom, 0}, {0, zoom}, {0, 0}};
            }
        }

        // update texture
        if (app->img->updated) {
            update_display_image(app->img);
            glBindTexture(GL_TEXTURE_2D, app->img->gl_txt);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, app->img->img.width(),
                app->img->img.height(), GL_RGBA, GL_FLOAT,
                app->img->display.data());
            app->img->updated = false;
        }

        // draw
        draw(win);

        // event hadling
        if (mouse_button || widgets_active) {
            glfwPollEvents();
        } else {
            glfwWaitEvents();
        }
    }
}

int main(int argc, char* argv[]) {
    // command line parameters
    std::vector<std::string> filenames = {"img.png"};
    float gamma = 2.2f;     // hdr gamma
    float exposure = 0.0f;  // hdr exposure
    bool filmic = false;    // filmic tone mapping
    bool diff = false;      // compute diffs
    bool lum_diff = false;  // compute luminance diffs
    bool quiet = false;     // quiet mode

    // command line params
    CLI::App parser("view images", "yimview");
    parser.add_option("--gamma,-g", gamma, "display gamma");
    parser.add_option("--exposure,-e", exposure, "display exposure");
    parser.add_flag("--filmic", filmic, "display filmic");
    parser.add_flag("--diff,-d", diff, "compute diff images");
    parser.add_flag("--luminance-diff,-D", lum_diff, "compute luminance diffs");
    parser.add_flag("--quiet,-q", quiet, "Print only errors messages");
    parser.add_option("images", filenames, "image filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // prepare application
    auto app = std::make_shared<app_state>();

    // loading images
    for (auto filename : filenames) {
        if (!quiet) std::cout << "loading " << filename << "\n";
        app->imgs.push_back(load_gimage(filename, exposure, gamma, filmic));
    }
    app->img = app->imgs.at(0);
    if (diff) {
        for (auto i = 0; i < app->imgs.size(); i++) {
            for (auto j = i + 1; j < app->imgs.size(); j++) {
                if (!quiet)
                    std::cout << "diffing " << app->imgs[i]->filename << " "
                              << app->imgs[j]->filename << "\n";
                auto diff = diff_gimage(app->imgs[i], app->imgs[j], !lum_diff);
                if (diff) app->imgs.push_back(diff);
            }
        }
    }

    // run ui
    run_ui(app);

    // done
    return 0;
}
