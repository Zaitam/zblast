cmake_minimum_required(VERSION 3.20)
project(zblast CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

add_compile_options(-Wall -Wextra)
add_compile_definitions("$<$<CONFIG:Debug>:_GLIBCXX_ASSERTIONS>")

add_executable(zblast main.cpp)
