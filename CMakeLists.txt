cmake_minimum_required(VERSION 2.8.4)
project(FASTQuick)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/bin/cfg.sh "#!/bin/bash\n")
file(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/bin/cfg.sh "CMAKE_CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"\n")

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick.sh
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/bin/cfg.sh ${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick_template.sh >${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick.sh
        COMMAND chmod 755 ${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick.sh
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/bin/cfg.sh")

add_custom_target(run ALL
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick.sh)

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/include/ /usr/local/include/)
if (NOT HTS_INCLUDE_DIRS)
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif ()

find_library(HTS_LIBRARIES hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if (NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif ()


find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

add_subdirectory(VerifyBamID)

set(STATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VerifyBamID/statgen/")

add_subdirectory(libbwa)
add_subdirectory(misc)
add_subdirectory(src)

#add_subdirectory(test)

enable_testing()
add_test(NAME Test1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/example
        COMMAND FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref ref.test.fa --out_prefix test_out_ref)
add_test(NAME Test2
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/example
        COMMAND FASTQuick align --fq_list fq.test.list --index_prefix test_out_ref --out_prefix test_out)

add_test(NAME Test3
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/example
        COMMAND bash -c "samtools sort -O BAM -o test_out.sorted.bam test_out.bam && samtools index test_out.sorted.bam && ${CMAKE_CURRENT_SOURCE_DIR}/bin/FASTQuick pop+con --DisableSanityCheck --BamFile test_out.sorted.bam --SVDPrefix ../resource/hapmap_3.3.b37.dat --Reference ref.test.fa --Output test_out")