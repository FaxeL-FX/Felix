list(APPEND FELIX_SOURCES
    Felix.cpp
    bmp/BMPWriter.cc
    math/math.cpp
    object/object.cpp
)

target_sources(Felix PRIVATE ${FELIX_SOURCES})
target_sources(FelixInfsim PRIVATE ${FELIX_SOURCES})