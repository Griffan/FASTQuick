AUTOMAKE_OPTIONS=foreign     

noinst_LIBRARIES=libvcf_statgen.a

libvcf_statgen_a_SOURCES= VcfFile.cpp        VcfFileWriter.cpp     VcfGenotypeFormat.cpp  VcfHeader.cpp  VcfRecord.cpp  VcfRecordField.cpp   VcfRecordGenotype.cpp  VcfSubsetSamples.cpp VcfFileReader.cpp  VcfGenotypeField.cpp  VcfGenotypeSample.cpp  VcfHelper.cpp  VcfRecordDiscardRules.cpp  VcfRecordFilter.cpp  VcfRecordInfo.cpp

libvcf_statgen_a_CXXFLAGS= -O3 -lz -I../general/   -Wall -std=c++0x   -mmmx -msse -msse2 -msse3 -pipe -Werror  -Wno-strict-overflow -I../include -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

