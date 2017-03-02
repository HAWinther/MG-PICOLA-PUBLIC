#!/bin/bash

fftinc="$HOME/local/include"
fftlib="$HOME/local/lib"
compiler="g++-mp-6"

$compiler -O3 main.cpp -I$fftinc -L$fftlib -lfftw3 -lfftw3_threads -fopenmp -o simplepofk
