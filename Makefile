TARGET = dsc
LIBS = -lgsl -lncurses -lm -lblas
CC = gcc
CFLAGS = -O0 -g -Wall -std=gnu11 -I/opt/local/include/ -fopenmp
LDFLAGS = -L/opt/local/lib/

OS := $(shell uname)
ifeq ($(OS), Darwin)
CC = gcc-mp-7
endif

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c cubature/hcubature.c libprogressbar/*.c inih/*.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o cubature/*.o libprogressbar/*.o inih/*.o
	-rm -f $(TARGET)
