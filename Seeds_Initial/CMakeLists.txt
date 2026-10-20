cmake_minimum_required(VERSION 3.12)
project(Seeds)
find_package(Imagine REQUIRED COMPONENTS Images)

add_executable(Seeds Seeds.cpp)
target_link_libraries(Seeds Imagine::Images)
