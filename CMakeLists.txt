cmake_minimum_required(VERSION 3.12)
project(Panorama)
list(APPEND CMAKE_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED COMPONENTS LinAlg Images)

add_executable(Panorama src/panorama.cpp)
target_link_libraries(Panorama Imagine::LinAlg Imagine::Images Imagine::Graphics)