list(APPEND FELIX_SOURCES
    Felix.cpp
    bmp/BMPWriter.cc
    math/math.cpp
    object/object.cpp
)

target_sources(Felix PRIVATE ${FELIX_SOURCES})
target_sources(FelixInfsim PRIVATE ${FELIX_SOURCES})

if(TARGET FelixMath)
    target_sources(FelixMath PRIVATE math/math.cpp)
    target_include_directories(FelixMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()
