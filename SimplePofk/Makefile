SHELL := /bin/bash

# Set compiler
CC = g++

# Options
# Choose among mass-assignment functions using -D_FUNCTION, where FUNCTION = NGP, CIC or TSC 
# If the input Gadget data is in kpc/h, use -D_UNITS, otherwise it is assumed that it is already in Mpc/h

OPTIONS = -D_SUBTRACTSHOTNOISE -D_NGP
#-D_UNITS

# Path to FFTW library
I = -I/home/topolski/include
L = -L/home/topolski/lib -lfftw3 -lfftw3_threads

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
