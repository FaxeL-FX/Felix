cmake_minimum_required(VERSION 3.20)

project(Felix)

add_executable(Felix)

add_executable(FelixInfsim)
target_compile_definitions(FelixInfsim PRIVATE infsimIsHere=1)

add_subdirectory(Felix)