cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fpic")

include_directories(../general/ ../)
set(libbam_a_SOURCES BamIndex.cpp GenericSamInterface.cpp PosList.cpp SamFileHeader.cpp SamHeaderPG.cpp SamHeaderSQ.cpp
        SamQuerySeqWithRefHelper.cpp SamRecordPool.cpp SamTags.cpp BamInterface.cpp PileupElementBaseQual.cpp
        SamCoordOutput.cpp SamFilter.cpp SamHeaderRecord.cpp SamHeaderTag.cpp SamRecord.cpp SamReferenceInfo.cpp
        SamValidation.cpp CigarHelper.cpp PileupElement.cpp SamFile.cpp SamHeaderHD.cpp SamHeaderRG.cpp SamInterface.cpp
        SamRecordHelper.cpp SamStatistics.cpp)

add_library(bam ${libbam_a_SOURCES})
target_link_libraries(bam general_statgen)

