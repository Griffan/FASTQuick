cmake_minimum_required(VERSION 2.8.4)

file(GLOB CXX_SOURCE_FILES *.cpp)

add_definitions(-DHAVE_PTHREAD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -Wall -fpic -pthread  -Wno-format-security -Wno-sign-compare -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

include_directories(../misc/vcf/ ${STATGEN_INCLUDE_DIR} ../misc/)

#add_library(FASTQuick STATIC ${CXX_SOURCE_FILES})
add_executable(FASTQuick ${CXX_SOURCE_FILES})

target_link_libraries(FASTQuick bwa mpu bam tabix vcf_statgen statgen misc VerifyBamID)
