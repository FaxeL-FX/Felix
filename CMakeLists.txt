cmake_minimum_required(VERSION 3.20)

project(Felix)

set(CMAKE_CXX_STANDARD 26)

add_library(msvc_utf8 INTERFACE)
if(WIN32)
    target_compile_options(msvc_utf8 INTERFACE "/utf-8")
endif()

add_library(Felix STATIC)
add_executable(Terminal)


target_link_libraries(Felix PRIVATE msvc_utf8)
target_link_libraries(Terminal PRIVATE msvc_utf8 Felix)

add_subdirectory(Felix)