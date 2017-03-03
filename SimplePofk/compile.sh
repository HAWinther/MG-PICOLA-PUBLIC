#!/bin/bash

fftinc="$HOME/local/include"
fftlib="$HOME/local/lib"

g++ -O3 -Wall main.cpp -I$fftinc -L$fftlib -lfftw3 -lfftw3_threads -fopenmp -o simplepofk

