SHELL := /bin/bash

# Set compiler
CC = mpicc-openmpi-gcc6

# Deconvolve the windowfunction?
OPTIONS = -D_DECONVOLVEWINDOWFUNCTION

I = -I$(HOME)/local/include
L = -L$(HOME)/local/lib -lm -lfftw3_mpi -lfftw3

C = -O3 -Wall $(OPTIONS)

TARGETS := simplepofk_mpi
all: $(TARGETS)

OBJS = main.c
main.o		: main.c io_gadget.h Makefile

simplepofk_mpi: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

