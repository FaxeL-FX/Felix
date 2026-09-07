list(APPEND FELIX_SOURCES
    Felix.cpp
    math/math.cpp
    object/object.cpp
    text/text_image.cpp
)

target_sources(Felix PRIVATE ${FELIX_SOURCES})
target_sources(Terminal PRIVATE Terminal.cpp ImageListBinary.cpp)

target_include_directories(Terminal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(TARGET FelixMath)
    target_sources(FelixMath PRIVATE math/math.cpp)
    target_include_directories(FelixMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(TARGET Bot)
    target_sources(Bot PRIVATE Bot.cpp ImageListBinary.cpp)
    target_link_libraries(Bot PRIVATE dpp Felix)
endif()
