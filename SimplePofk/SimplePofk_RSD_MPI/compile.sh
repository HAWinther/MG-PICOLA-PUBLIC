#!/bin/bash
mpicc-openmpi-gcc6 -O3 main.c -I$HOME/local/include -L$HOME/local/lib -lfftw3 -lfftw3_threads -fopenmp -o simplepofk_mpi
mpicc-openmpi-gcc6 -O3 rsd_main.c -I$HOME/local/include -L$HOME/local/lib -lfftw3 -lfftw3_threads -fopenmp -o simplepofk_rsd_mpi
