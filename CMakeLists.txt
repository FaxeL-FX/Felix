cmake_minimum_required(VERSION 3.20)

project(Felix)

set(CMAKE_CXX_STANDARD 26)

add_library(msvc_utf8 INTERFACE)
if(WIN32)
    target_compile_options(msvc_utf8 INTERFACE "/utf-8")
endif()

add_executable(Felix)


target_link_libraries(Felix PRIVATE msvc_utf8)

add_subdirectory(Felix)