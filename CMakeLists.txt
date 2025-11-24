cmake_minimum_required(VERSION 3.20)

project(Felix)

add_library(msvc_utf8 INTERFACE)
if(WIN32)
    target_compile_options(msvc_utf8 INTERFACE "/utf-8")
endif()

add_executable(Felix)

add_executable(FelixInfsim)
target_compile_definitions(FelixInfsim PRIVATE infsimIsHere=1)

target_link_libraries(Felix PRIVATE msvc_utf8)
target_link_libraries(FelixInfsim PRIVATE msvc_utf8)

add_subdirectory(Felix)