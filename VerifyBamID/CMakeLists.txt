cmake_minimum_required(VERSION 2.8.4)
project(VerifyBamID)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-variable -Wno-unused-result -Wno-format-security -Wno-overloaded-virtual")
add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__ZLIB_AVAILABLE__)


find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()
find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(SSLLIB ssl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT SSLLIB)
    message(FATAL_ERROR "libssl library not found")
endif()
find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()


add_subdirectory(statgen)

set(LIBSTATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/statgen")

add_subdirectory(samtools)

add_subdirectory(libVcf)

set(SOURCE_FILES vb2Main.cpp SVDcalculator.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp params.cpp )#SVDcalculator.cpp SVDcalculator.h)

include_directories(statgen ${HTS_INCLUDE_DIRS} samtools libVcf Eigen)

add_library(VerifyBamID ${SOURCE_FILES})

target_link_libraries(VerifyBamID statgen Vcf ${HTS_LIBRARIES} samtools ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB} ${SSLLIB} ${CRYPTOLIB})

enable_testing()

add_test( NAME myTest1
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMAND VerifyBamID --DisableSanityCheck --UDPath resource/test/hapmap_3.3.b37.dat.UD --BamFile resource/test/test.bam --BedPath resource/test/hapmap_3.3.b37.dat.bed --MeanPath resource/test/hapmap_3.3.b37.dat.mu --Reference resource/test/chr20.fa.gz )
#add_test( NAME myTest2
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          COMMAND sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run.plot.sh -i ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap_3.3.b37.dat.V -o ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap -r 1000g -g grey)
