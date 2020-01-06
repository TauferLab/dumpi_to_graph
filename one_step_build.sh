#!/usr/bin/env bash

export CC=mpicc
export CXX=mpicxx
rm -rf $(pwd)/build && mkdir -p $(pwd)/build && cd build
cmake -DDUMPI_ROOT=$HOME/repos/sst-dumpi/build ..
make -j
