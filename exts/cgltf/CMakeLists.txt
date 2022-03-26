cmake_minimum_required(VERSION 3.12.0)

project(cgltf)

add_library(cgltf STATIC 
  cgltf/cgltf.c
  cgltf/cgltf.h
  cgltf/cgltf_write.c
  cgltf/cgltf_write.h
)

target_include_directories(cgltf PUBLIC . cgltf)

add_library(cgltf::cgltf ALIAS cgltf)
