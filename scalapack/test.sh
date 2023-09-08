#! /bin/bash

mkdir build
cd build
cmake ..
make -j
cd ../scalapack-driver/
make
mpirun -n 64 ./bin/pdqrdriver
cat QR.out
