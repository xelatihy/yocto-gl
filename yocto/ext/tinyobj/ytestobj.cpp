#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_utils.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "../yocto/ext/tinyobj/tiny_obj_loader.h"
#include <chrono>

int main(int argc, char** argv) {
    auto filename = std::string(argv[1]);
    auto dirname = ygl::path_dirname(filename);
    auto clock = std::chrono::high_resolution_clock();

    auto start_yobj = clock.now();
    auto obj = ygl::load_obj(filename, true);
    auto end_yobj = clock.now();  

    auto start_tobj = clock.now();
    auto attrib = tinyobj::attrib_t();
    auto shapes = std::vector<tinyobj::shape_t>();
    auto materials = std::vector<tinyobj::material_t>();
    auto err = std::string();
    if(!tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str(), dirname.c_str())) throw std::runtime_error("could not load obj");
    auto end_tobj = clock.now();

    auto elapsed_tobj = std::chrono::duration<double, std::milli>(end_tobj - start_tobj).count();
    auto elapsed_yobj = std::chrono::duration<double, std::milli>(end_yobj - start_yobj).count();

    printf("tinyobj:   %lf %d\n", elapsed_tobj, (int)shapes.size());
    printf("yocto_obj: %lf %d\n", elapsed_yobj, (int)obj->objects.size());
}