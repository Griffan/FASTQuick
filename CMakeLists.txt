cmake_minimum_required(VERSION 2.8.4)
project(FASTQuick)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3  -std=c++0x -Wno-unused-variable -Wno-unused-result -Wno-format-security")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

find_library(HTS_LIBRARIES hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

add_subdirectory(VerifyBamID)

set(STATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VerifyBamID/statgen/")

add_subdirectory(libbwa)
add_subdirectory(libmpu)
add_subdirectory(misc)
add_subdirectory(src)
#add_subdirectory(test)