#!/bin/bash

cd ..
PREFIX=$(pwd)/$(uname)/$(uname -m)
LIB=$PREFIX/lib
rm -rf LIB/*
cd build

rm -rf Ipopt Thirparty build

curl -o coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod a+x coinbrew
export MACOSX_DEPLOYMENT_TARGET=$MIN_MACOSX_VERSION
./coinbrew fetch Ipopt@3.14 --without-asl --without-hsl --with-mumps --with-metis 
./coinbrew build Ipopt@3.14 --prefix $PREFIX --disable-openmp --without-asl --without-hsl --with-mumps --with-metis --enable-static=yes --enable-shared=no --parallel-jobs 4 
