SHELL := /bin/bash

# Set compiler
CC = g++-mp-6

# Options
OPTIONS = -D_SUBTRACTSHOTNOISE -D_CIC

# Path to FFTW library
I = -I$(HOME)/local/include
L = -L$(HOME)/local/lib -lfftw3 -lfftw3_threads

C = -O3 -Wall -fopenmp $(OPTIONS)

TARGETS := simplepofk
all: $(TARGETS)

OBJS = main.cpp
main.o: io_ascii.h  io_gadget.h io_ramses.h Makefile

simplepofk: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o
