cmake_minimum_required(VERSION 3.12.0)

project(tinyexr)

add_library(tinyexr STATIC 
  tinyexr/tinyexr.cpp
)

target_include_directories(tinyexr PUBLIC . tinyexr)

add_library(tinyexr::tinyexr ALIAS tinyexr)
