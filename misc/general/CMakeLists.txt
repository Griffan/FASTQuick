cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fpic")
add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__ZLIB_AVAILABLE__)
add_definitions(-D_FILE_OFFSET_BITS=64)
include_directories(../)
set(libgeneral_a_SOURCES BaseAsciiMap.cpp   BgzfFileTypeRecovery.cpp  CigarRoller.cpp   FileType.cpp
        GzipFileType.cpp  IndexBase.cpp  MathMatrix.cpp      MemoryMap.cpp  StatGenStatus.cpp  StringBasics.cpp
        UncompressedFileType.cpp BaseUtilities.cpp  CharBuffer.cpp            Error.cpp         Generic.cpp
        GzipHeader.cpp    InputFile.cpp  MathVector.cpp      PhoneHome.cpp  STLUtilities.cpp   StringHash.cpp
        BgzfFileType.cpp   Cigar.cpp                 ErrorHandler.cpp  GenomeSequence.cpp  Hash.cpp
        IntArray.cpp   MemoryMapArray.cpp  Sort.cpp       StringArray.cpp    Tabix.cpp)

add_library(general_statgen ${libgeneral_a_SOURCES})
target_link_libraries(general_statgen z)

