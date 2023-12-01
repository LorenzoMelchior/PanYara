cat CMakeLists.txt
cmake_minimum_required (VERSION 3.14)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wpedantic -fPIC")

project(PanYara)

find_package(SeqAn3 QUIET REQUIRED HINTS /include/seqan3/build_system)
find_package(IVio QUIET REQUIRED HINTS /include/IVio)


add_executable(pan_yara_indexer pan_indexer.cpp
                                pibf.cpp)

add_executable(pan_yara_mapper pan_mapper.cpp
                               pibf.cpp)

target_link_libraries(pan_yara_indexer PUBLIC ivio::ivio seqan3::seqan3)
target_link_libraries(pan_yara_mapper PUBLIC ivio::ivio seqan3::seqan3)
