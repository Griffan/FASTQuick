AUTOMAKE_OPTIONS=foreign     

noinst_LIBRARIES=libmisc.a

#SUBDIRS=general vcf bam

libmisc_a_SOURCES= faidx.c params.cpp razf.c knetfile.c bgzf.c

#libmisc_a_LIBADD= ./general/libgeneral.a ./vcf/libvcf_statgen.a ./bam/libbam.a

libmisc_a_CXXFLAGS= -O3  -I../misc/vcf/ -I../misc/general/ -lz  -Wall -std=c++0x  -mmmx -msse -msse2 -msse3 -pipe  -Wno-strict-overflow -I../include -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

libmisc_a_CFLAGS=  -O3 -lz  -Wall  -Wno-unused-result -mmmx -msse -msse2 -msse3 -pipe  -I../include -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

SUBDIRS=general vcf bam
