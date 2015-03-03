AUTOMAKE_OPTIONS=foreign     

bin_PROGRAMS=FASTQuick

FASTQuick_SOURCES=FASTQuick.cpp BwtIndexer.cpp BwtMapper.cpp RefBuilder.cpp StatCollector.cpp

FASTQuick_LDFLAGS= 

FASTQuick_LDADD= ../libmpu/libtabix.a ../libmpu/libmpu.a ../libbwa/libbwa.a ../misc/libmisc.a ../libStatGen/libStatGen.a

FASTQuick_CXXFLAGS= -O4   -Wall -std=c++11  -DHAVE_PTHREAD -fpermissive -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -mmmx -msse -msse2 -msse3 -I../libStatGen/include/ -L../libStatGen/ -L../libmpu/ -I../misc/  -lStatGen -lmpu 

