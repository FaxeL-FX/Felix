cmake_minimum_required(VERSION 3.20)

project(Felix)

set(CMAKE_CXX_STANDARD 26)

#find_package(dpp)
#set(BUILD_SHARED_LIBS OFF)
#add_subdirectory("D:\\D\\Projects\\VS\\DPP-10.1.6" DPP_BUILD_DIR)
set(FELIX_DPP_LIBRARY_PATH "" CACHE STRING "Path to discord api library [DPP-10.1.6]")
if(FELIX_DPP_LIBRARY_PATH)
    add_subdirectory(${FELIX_DPP_LIBRARY_PATH} DPP_BUILD_DIR)
else()
    message("Felix discord bot is disabled")
endif()

add_library(msvc_utf8 INTERFACE)
if(WIN32)
    target_compile_options(msvc_utf8 INTERFACE "/utf-8")
endif()

add_library(Felix STATIC)
add_executable(Terminal)
if(FELIX_DPP_LIBRARY_PATH)
    add_executable(Bot)
endif()

target_link_libraries(Felix PRIVATE msvc_utf8)
target_link_libraries(Terminal PRIVATE msvc_utf8 Felix)

add_subdirectory(Felix)
