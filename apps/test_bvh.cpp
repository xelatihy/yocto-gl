//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

// #define YGL_USESTL

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#define NANORT_IMPLEMENTATION
#define NANORT_LOG
#include "ext/nanort.h"

#define YGL_BVH_LOG_RAYS
#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_obj.h"

void collect_boxes(yb_scene* scene, int shape_id, int nodeid, int max_depth,
                   int depth, ym_vector<ym_range3f>& bbox) {
    yb__bvh* bvh = (shape_id >= 0) ? scene->shapes[shape_id].bvh : scene->bvh;
    if (depth > max_depth) return;
    yb__bvhn* node = &bvh->nodes[nodeid];
    bbox += node->bbox;
    if (!node->isleaf) {
        for (int i = 0; i < node->count; i++) {
            collect_boxes(scene, shape_id, node->start + i, max_depth,
                          depth + 1, bbox);
        }
    } else {
        if (shape_id < 0) {
            for (int i = 0; i < node->count; i++) {
                collect_boxes(scene, bvh->sorted_prim[node->start + i], 0,
                              max_depth, depth + 1, bbox);
            }
        }
    }
}

void collect_boxes(nanort::BVHAccel* bvh, int nodeid, int max_depth, int depth,
                   ym_vector<ym_range3f>& bbox) {
    if (depth > max_depth) return;
    const nanort::BVHNode* node = &bvh->GetNodes()[nodeid];
    bbox += {ym_vec3f(node->bmin), ym_vec3f(node->bmax)};
    if (!node->flag) {
        collect_boxes(bvh, node->data[0], max_depth, depth + 1, bbox);
        collect_boxes(bvh, node->data[1], max_depth, depth + 1, bbox);
    }
}

ym_vector<ym_range3f> collect_boxes(yb_scene* scene, int max_depth) {
    ym_vector<ym_range3f> boxes;
    collect_boxes(scene, -1, 0, max_depth, 0, boxes);
    return boxes;
}

ym_vector<ym_range3f> collect_boxes(nanort::BVHAccel* bvh, int max_depth) {
    ym_vector<ym_range3f> boxes;
    collect_boxes(bvh, 0, max_depth, 0, boxes);
    return boxes;
}

void save_bvh_obj(const char* filename, const ym_vector<ym_range3f>& bboxes) {
    yo_scene* scene = new yo_scene();
    for (auto&& bbox : bboxes) {
        scene->shapes.push_back(yo_shape());
        yo_shape* shape = &scene->shapes.back();
        shape->etype = 2;
        shape->nverts = 8;
        shape->pos = {
            {bbox.min.x, bbox.min.y, bbox.min.z},
            {bbox.min.x, bbox.min.y, bbox.max.z},
            {bbox.min.x, bbox.max.y, bbox.min.z},
            {bbox.min.x, bbox.max.y, bbox.max.z},
            {bbox.max.x, bbox.min.y, bbox.min.z},
            {bbox.max.x, bbox.min.y, bbox.max.z},
            {bbox.max.x, bbox.max.y, bbox.min.z},
            {bbox.max.x, bbox.max.y, bbox.max.z},
        };
        shape->nelems = 12;
        shape->elem = {0, 1, 1, 3, 3, 2, 2, 0, 4, 5, 5, 7,
                       7, 6, 6, 4, 0, 4, 1, 5, 2, 6, 3, 7};
    }
    yo_save_obj(filename, scene, true);
    delete scene;
}

void convert_bvh_node(nanort::BVHAccel* nbvh, int ni, yb__bvh* bvh, int i,
                      int* nnodes) {
    const nanort::BVHNode* nn = &nbvh->GetNodes()[ni];
    yb__bvhn* n = &bvh->nodes[i];
    n->axis = nn->axis;
    n->bbox = ym_range3f{ym_vec3f{nn->bmin}, ym_vec3f{nn->bmax}};
    n->count = 2;
    n->isleaf = nn->flag;
    if (n->isleaf) {
        n->start = nn->data[1];
        n->count = nn->data[0];
    } else {
        n->start = *nnodes;
        n->count = 2;
        *nnodes += 2;
        convert_bvh_node(nbvh, nn->data[0], bvh, n->start + 0, nnodes);
        convert_bvh_node(nbvh, nn->data[1], bvh, n->start + 1, nnodes);
    }
}

void print_bvh(nanort::BVHAccel* bvh, int nid, int d) {
    static char axes[3] = {'x', 'y', 'z'};
    const nanort::BVHNode* n = &bvh->GetNodes()[nid];
    for (int i = 0; i < d; i++) printf(" ");
    printf("%4d", nid);
    if (n->flag) {
        printf(" / ");
        for (int i = 0; i < n->data[0]; i++) printf("%d ", n->data[1] + i);
    } else
        printf(" / %c ", axes[n->axis]);
    printf("\n");

    if (!n->flag) {
        print_bvh(bvh, n->data[0], d + 1);
        print_bvh(bvh, n->data[1], d + 1);
    }
}

void print_bvh(yb_scene* scene, int sid, int nid, int d) {
    yb__bvh* bvh = (sid >= 0) ? scene->shapes[sid].bvh : scene->bvh;
    static char axes[3] = {'x', 'y', 'z'};
    const yb__bvhn* n = &bvh->nodes[nid];
    for (int i = 0; i < d; i++) printf(" ");
    printf("%4d", nid);
    if (n->isleaf) {
        if (sid >= 0) {
            printf(" / ");
            for (int i = 0; i < n->count; i++) printf(" %d", n->start + i);
        } else {
            for (int i = 0; i < n->count; i++)
                print_bvh(scene, bvh->sorted_prim[n->start + i], 0, d + 1);
        }
    } else {
        printf(" / %c ", axes[n->axis]);
    }
    printf("\n");

    if (!n->isleaf) {
        for (int i = 0; i < n->count; i++)
            print_bvh(scene, sid, n->start + i, d + 1);
    }
}

int main(int argc, const char** argv) {
    // command line
    yc_parser* parser = yc_init_parser(argc, argv, "test obj");
    bool flatten = yc_parse_optb(parser, "--flatten", "-f",
                                 "flatten scene to one shape", false);
    const char* filename =
        yc_parse_args(parser, "scene", "scene filename", 0, true);
    yc_done_parser(parser);

    // paths
    char dirname[512], basename[512];
    yc_split_path(filename, dirname, basename, 0);

    // scene and camera
    ym_timer yocto_load_timer = ym_timer();
    yo_scene* scene =
        yo_load_obj(filename, true, true);  // make camera if not there
    if (!scene->cameras.size()) {
        // find scene bounds
        ym_range3f bbox = ym_invalid_range3f;
        for (int i = 0; i < scene->shapes.size(); i++) {
            yo_shape* shape = &scene->shapes[i];
            for (int j = 0; j < shape->nverts; j++) bbox += shape->pos[j];
        }
        ym_vec3f bbox_center = ym_rcenter(bbox);
        ym_vec3f bbox_size = ym_rsize(bbox);
        float bbox_msize = fmax(bbox_size.x, fmax(bbox_size.y, bbox_size.z));
        // create camera
        scene->cameras.resize(1);
        yo_camera* cam = &scene->cameras[0];
        // set up camera
        ym_vec3f camera_dir = {1, 0.4f, 1};
        cam->from = camera_dir * bbox_msize + bbox_center;
        cam->to = bbox_center;
        cam->up = ym_vec3f{0, 1, 0};
        cam->width = 16.0f / 9.0f;
        cam->height = 1;
        cam->aperture = 0;
    }
    ym_frame3f camera_xf = ym_lookat_xform3(
        scene->cameras[0].from, scene->cameras[0].to, scene->cameras[0].up);
    double yocto_load_elapsed = yocto_load_timer.elapsed();
    printf("load obj: %d shapes in %lg ms\n", (int)scene->shapes.size(),
           yocto_load_elapsed);

    // flatten if needed
    if (flatten) {
        yo_shape flattened;
        flattened.name = "flattened";
        flattened.matname = "";
        flattened.matid = -1;
        flattened.etype = 3;
        for (int i = 0; i < scene->shapes.size(); i++) {
            int offset = (int)flattened.pos.size();
            yo_shape* shape = &scene->shapes[i];
            if (shape->etype != 3) continue;
            for (int j = 0; j < shape->pos.size(); j++)
                flattened.pos += shape->pos[j];
            for (int j = 0; j < shape->elem.size() / 3; j++) {
                flattened.elem += {shape->elem[j * 3 + 0] + offset,
                                   shape->elem[j * 3 + 1] + offset,
                                   shape->elem[j * 3 + 2] + offset};
            }
        }
        flattened.nelems = (int)flattened.elem.size() / 3;
        flattened.nverts = (int)flattened.pos.size();
        scene->shapes.clear();
        scene->shapes += flattened;
    }

    // resolution
    int w = 1280, h = 720, ns = 1, nns = 1;
    float cw = 1.77777777f, ch = 1;

    // build bvh
    ym_timer yocto_bvh_timer = ym_timer();
    yb_scene* yocto_bvh = yb_init_scene((int)scene->shapes.size());
    for (int i = 0; i < scene->shapes.size(); i++) {
        yo_shape* shape = &scene->shapes[i];
        yb_set_shape(yocto_bvh, i, shape->xform, shape->nelems,
                     shape->elem.data(), shape->etype, shape->nverts,
                     shape->pos.data(), shape->radius.data());
    }
    yb_build_bvh(yocto_bvh, 0);
    double yocto_bvh_elapsed = yocto_bvh_timer.elapsed();
    printf("build yocto: %d shapes in %lg ms\n", (int)scene->shapes.size(),
           yocto_bvh_elapsed);

    // nanort --- build bvh
    ym_timer nanort_bvh_timer = ym_timer();
    nanort::BVHBuildOptions nanort_options;
    nanort::BVHAccel nanort_bvh;
    ym_vector<ym_vec3f> nanort_pos;
    ym_vector<ym_vec3i> nanort_elem;
    for (int i = 0; i < scene->shapes.size(); i++) {
        int offset = (int)nanort_pos.size();
        yo_shape* shape = &scene->shapes[i];
        if (shape->etype != 3) continue;
        for (int j = 0; j < shape->pos.size(); j++) nanort_pos += shape->pos[j];
        for (int j = 0; j < shape->elem.size() / 3; j++) {
            nanort_elem += {shape->elem[j * 3 + 0] + offset,
                            shape->elem[j * 3 + 1] + offset,
                            shape->elem[j * 3 + 2] + offset};
        }
    }
    nanort_bvh.Build(&nanort_pos.data()->x,
                     (unsigned int*)&nanort_elem.data()->x,
                     (int)nanort_elem.size(), nanort_options);
    double nanort_bvh_elapsed = nanort_bvh_timer.elapsed();
    printf("nanort build bvh: %d shapes in %lg ms\n", 0, nanort_bvh_elapsed);

    // copy of the nanort bvh
    //    yb_scene_bvh* cbvh = convert_bvh(&nanort_bvh,
    //    (int*)nanort_elem.data(), (float*)nanort_pos.data());
    //    print_bvh(&nanort_bvh, 0, 0);
    //    print_bvh(cbvh->shapes[0], 0, 0);

    // save obj
    save_bvh_obj(
        (ym_string(dirname) + ym_string(basename) + ".yocto_bvh.obj").c_str(),
        collect_boxes(yocto_bvh, 8));
    save_bvh_obj(
        (ym_string(dirname) + ym_string(basename) + ".nanort_bvh.obj").c_str(),
        collect_boxes(&nanort_bvh, 7));

    // raycast
    ym_timer yocto_render_timer = ym_timer();
    ym_vector<ym_vec3f> img(w * h);
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            ym_vec3f c = ym_zero3f;
            for (int s = 0; s < ns; s++) {
                ym_vec2f uv = {(i + 0.5f) / w, (h - j + 0.5f) / h};
                ym_vec3f q = {cw * (uv.x - 0.5f), ch * (uv.y - 0.5f), -1};
                ym_ray3f ray = {
                    ym_transform_point(camera_xf, ym_zero3f),
                    ym_transform_direction(camera_xf, ym_normalize(q))};
                int sid, eid;
                float ray_t;
                ym_vec2f euv;
                if (yb_intersect_first(yocto_bvh, ray, &ray_t, &sid, &eid,
                                       &euv)) {
                    int* f = scene->shapes[sid].elem.data() + 3 * eid;
                    ym_vec3f* pos = scene->shapes[sid].pos.data();
                    ym_vec3f norm = ym_normalize(
                        ym_cross(pos[f[1]] - pos[f[0]], pos[f[2]] - pos[f[0]]));
                    c += ym_vec3f(-ym_dot(ray.d, norm));
                }
            }
            img[j * w + i] = c / ns;
        }
    }
    double yocto_render_elapsed = yocto_render_timer.elapsed();
    printf("render image: %dx%d@%d in %lg ms\n", w, h, ns,
           yocto_render_elapsed);

    // nanort --- raycast
    ym_timer nanort_render_timer = ym_timer();
    ym_vector<ym_vec3f> nanort_img(w * h);
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            ym_vec3f c = ym_zero3f;
            for (int s = 0; s < ns; s++) {
                ym_vec2f uv = {(i + 0.5f) / w, (h - j + 0.5f) / h};
                ym_vec3f q = {cw * (uv.x - 0.5f), ch * (uv.y - 0.5f), -1};
                ym_ray3f ray = {
                    ym_transform_point(camera_xf, ym_zero3f),
                    ym_transform_direction(camera_xf, ym_normalize(q))};

                nanort::Ray nanort_ray = nanort::Ray();
                nanort_ray.minT = ray.tmin;
                nanort_ray.maxT = ray.tmax;
                nanort_ray.org[0] = ray.o.x;
                nanort_ray.org[1] = ray.o.y;
                nanort_ray.org[2] = ray.o.z;
                nanort_ray.dir[0] = ray.d.x;
                nanort_ray.dir[1] = ray.d.y;
                nanort_ray.dir[2] = ray.d.z;

                nanort::Intersection isect;
                isect.t = 1.0e+30f;

                nanort::BVHTraceOptions traceOptions;
                if (nanort_bvh.Traverse(&isect, (float*)nanort_pos.data(),
                                        (unsigned int*)nanort_elem.data(),
                                        nanort_ray, traceOptions)) {
                    int* f = (int*)(nanort_elem.data() + isect.faceID);
                    ym_vec3f* pos = nanort_pos.data();
                    ym_vec3f norm = ym_normalize(
                        ym_cross(pos[f[1]] - pos[f[0]], pos[f[2]] - pos[f[0]]));
                    c += ym_vec3f(-ym_dot(ray.d, norm));
                }
            }
            nanort_img[j * w + i] = c / ns;
        }
    }
    double nanort_render_elapsed = nanort_render_timer.elapsed();
    printf("render image: %dx%d@%d in %lg ms\n", w, h, ns,
           nanort_render_elapsed);

    // stats
    printf("yocto bvh stats\n");
    int nprims, nleaves, ninternals, min_depth, max_depth;
    yb_compute_bvh_stats(yocto_bvh, false, &nprims, &ninternals, &nleaves,
                         &min_depth, &max_depth);
    printf("scene bvh nodes: %d\n", nleaves + ninternals);
    printf("scene bvh prims: %d\n", nprims);
    printf("scene bvh leaves: %d\n", nleaves);
    printf("scene bvh internals: %d\n", ninternals);
    printf("scene bvh depth: %d - %d\n", min_depth, max_depth);
    yb_compute_bvh_stats(yocto_bvh, true, &nprims, &ninternals, &nleaves,
                         &min_depth, &max_depth);
    printf("full bvh nodes: %d\n", nleaves + ninternals);
    printf("full bvh prims: %d\n", nprims);
    printf("full bvh leaves: %d\n", nleaves);
    printf("full bvh internals: %d\n", ninternals);
    printf("full bvh depth: %d - %d\n", min_depth, max_depth);
    printf("\nnanort bvh stats\n");
    printf("nanort bvh nodes: %d\n", (int)nanort_bvh.GetNodes().size());
    printf("nanort bvh prims: %d\n", (int)nanort_bvh.GetIndices().size());
    printf("nanort bvh leaves: %d\n", nanort_bvh.GetStatistics().numLeafNodes);
    printf("nanort bvh internals: %d\n",
           nanort_bvh.GetStatistics().numBranchNodes);
    printf("nanort bvh depth: %d\n", nanort_bvh.GetStatistics().maxTreeDepth);

    // logs
    int nrays, nbbox_inters, npoint_inters, nline_inters, ntriangle_inters;
    yb_get_ray_log(&nrays, &nbbox_inters, &npoint_inters, &nline_inters,
                   &ntriangle_inters);
    printf("\nray interssection log\n");
    printf("nrays: %d\n", nrays);
    printf("nbbox: %d\n", nbbox_inters);
    printf("npoints: %d\n", npoint_inters);
    printf("nlines: %d\n", nline_inters);
    printf("ntriangles: %d\n", ntriangle_inters);
    printf("\nanort ray interssection log\n");
    printf("nrays: %d\n", nanort::log_nrays);
    printf("nbbox: %d\n", nanort::log_nbbox_inters);
    printf("npoints: %d\n", nanort::log_npoint_inters);
    printf("nlines: %d\n", nanort::log_nline_inters);
    printf("ntriangles: %d\n", nanort::log_ntriangle_inters);

    // save
    ym_string yocto_img_filename =
        ym_string(dirname) + ym_string(basename) + ".yocto.hdr";
    stbi_write_hdr(yocto_img_filename.c_str(), w, h, 3, &img.data()->x);
    ym_string nanort_img_filename =
        ym_string(dirname) + ym_string(basename) + ".nanort.hdr";
    stbi_write_hdr(nanort_img_filename.c_str(), w, h, 3, &nanort_img.data()->x);
}
