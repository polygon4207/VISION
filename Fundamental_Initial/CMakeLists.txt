cmake_minimum_required(VERSION 3.12)
project(Fundamental)
set(CMAKE_CXX_STANDARD 20)

find_package(Imagine REQUIRED COMPONENTS LinAlg Images)

add_executable(Fundamental
        Fundamental.cpp
        Imagine/SIFT_VL.cpp
        Imagine/vl/generic.c Imagine/vl/host.c Imagine/vl/imop.c Imagine/vl/sift.c)
target_link_libraries(Fundamental Imagine::LinAlg Imagine::Images)
