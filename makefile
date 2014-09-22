TARGET = FastqA
LIBS = -lz -lStatGen -L/net/wonderland/home/fanzhang/WorkingSpace/lib -lpthread 
INCLUDE	= -IlibStatGen-master/include/ 
CC = g++
CFLAGS = -O3   -lStatGen -g   -Wall -std=c++11  -DHAVE_PTHREAD -fpermissive

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.cpp libbwa/*.c misc/*.c))
HEADERS = $(wildcard *.h  libbwa/*.h misc/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE)   -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -O3  -g  -Wall -std=c++11 -DHAVE_PTHREAD -mmmx -msse -msse2 -msse3   -fpermissive $(LIBS) $(INCLUDE) -o $@

clean:
	-rm -f *.o ./libbwa/*.o 
	-rm -f $(TARGET)
