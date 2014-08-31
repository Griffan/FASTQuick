TARGET = FastqA
LIBS = -lz -lStatGen -L/net/wonderland/home/fanzhang/WorkingSpace/lib -lpthread
INCLUDE	= -IlibStatGen-master/general/
CC = g++
CFLAGS = -O3   -g   -Wall -std=c++11  -DHAVE_PTHREAD -fpermissive

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.cpp samtools-faidx/*.c libbwa/*.c ))
HEADERS = $(wildcard *.h samtools-faidx/*.h libbwa/*.h )

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE)  -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -O3  -g -std=c++11 -Wall -DHAVE_PTHREAD -mmmx -msse -msse2 -msse3 -fpermissive $(LIBS) $(INCLUDE) -o $@

clean:
	-rm -f *.o ./libbwa/*.o ./samtools-faidx/*.o
	-rm -f $(TARGET)
