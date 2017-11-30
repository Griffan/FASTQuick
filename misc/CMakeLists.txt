cmake_minimum_required(VERSION 2.8.4)
add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__ZLIB_AVAILABLE__)
add_definitions(-D_FILE_OFFSET_BITS=64)
add_subdirectory(general)
add_subdirectory(vcf)
add_subdirectory(bam)

include_directories(general)
set(libmisc_a_SOURCES faidx.c params.cpp razf.c knetfile.c bgzf.c)
add_library(misc ${libmisc_a_SOURCES})

