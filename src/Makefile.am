AUTOMAKE_OPTIONS=foreign     

bin_PROGRAMS=FASTQuick

FASTQuick_SOURCES=FASTQuick.cpp BwtIndexer.cpp BwtMapper.cpp RefBuilder.cpp StatCollector.cpp InsertSizeEstimator.cpp PopulationIdentifier.cpp ContaminationEstimator.cpp

FASTQuick_LDFLAGS= 

FASTQuick_LDADD= ../libbwa/libbwa.a ../misc/libmisc.a ../misc/general/libgeneral.a ../libmpu/libmpu.a ../misc/bam/libbam.a  ../misc/libmisc.a ../libmpu/libtabix.a ../misc/general/libgeneral.a ../misc/vcf/libvcf_statgen.a ../misc/general/libgeneral.a 

FASTQuick_CXXFLAGS=  -O3 -g -std=c++0x  -DHAVE_PTHREAD  -mmmx -msse -msse2 -msse3  -I../misc/vcf/ -I../misc/general/ -I../misc/ -lz -lpthread -pipe  -Wno-strict-overflow -I../include -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS -Wno-unused-variable -Wno-unused-function
