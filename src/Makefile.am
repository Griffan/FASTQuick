AUTOMAKE_OPTIONS=foreign     

bin_PROGRAMS=FastqA

FastqA_SOURCES=FastqA.cpp BwtIndexer.cpp BwtMapper.cpp RefBuilder.cpp StatCollector.cpp

FastqA_LDFLAGS= 

FastqA_LDADD= ../libmpu/libtabix.a ../libmpu/libmpu.a ../libbwa/libbwa.a ../misc/libmisc.a -lStatGen -lmpu -ltabix

FastqA_CXXFLAGS= -O3   -g   -Wall -std=c++11  -DHAVE_PTHREAD -fpermissive -Wno-unused-but-set-variable -D_GLIBCXX_DEBUG -Wno-unused-variable -mmmx -msse -msse2 -msse3 -I../libStatGen/include/ -L../libStatGen/ -L../libmpu/ -I../misc/

