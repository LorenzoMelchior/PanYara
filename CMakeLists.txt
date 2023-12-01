cmake_minimum_required (VERSION 3.14)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PanYara)

include_directories(${PROJECT_SOURCE_DIR}/include)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan3/include/seqan3/version.hpp")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan3/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan3/build_system
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
       ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan3/build_system
       ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn3 library provided with the PAN-YARA source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

find_package(SeqAn3 QUIET REQUIRED CONFIG)
find_package(IVio PATHS include/IVio)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_executable(pan_yara_indexer src/pan_indexer.cpp
                                src/pibf.cpp)

add_executable(pan_yara_mapper src/pan_mapper.cpp
                               src/pibf.cpp)

target_link_libraries(pan_yara_indexer PUBLIC ivio::ivio seqan3::seqan3)
target_link_libraries(pan_yara_mapper PUBLIC ivio::ivio seqan3::seqan3)
