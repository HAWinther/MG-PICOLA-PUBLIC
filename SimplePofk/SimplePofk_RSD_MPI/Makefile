SHELL := /bin/bash

# Set compiler
CC = mpicc-openmpi-gcc6 -std=c99

# Deconvolve the windowfunction?
OPTIONS = -D_DECONVOLVEWINDOWFUNCTION -D_READVELOCITY

I = -I$(HOME)/local/include
L = -L$(HOME)/local/lib -lfftw3_mpi -lfftw3

C = -O3 -Wall $(OPTIONS)

TARGETS := simplepofk_rsd_mpi
all: $(TARGETS)

# OBJECT FILES
OBJS = rsd_main.c

# DEPENDENCIES
main.o		: rsd_main.c io_gadget.h Makefile

# HEADERS
HEADERS = io_gadget.h

simplepofk_rsd_mpi: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

