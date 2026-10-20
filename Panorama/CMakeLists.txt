cmake_minimum_required(VERSION 3.12)
project(Panorama)
find_package(Imagine REQUIRED COMPONENTS LinAlg Images)

add_executable(Panorama Panorama.cpp)
target_link_libraries(Panorama Imagine::LinAlg Imagine::Images)
