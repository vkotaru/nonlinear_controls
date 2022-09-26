#!/bin/bash

mkdir tmp
cd tmp

git clone --recursive https://github.com/osqp/osqp
git clone --recursive https://github.com/embotech/ecos.git
git clone --recursive https://github.com/robotology/osqp-eigen.git

cd osqp
mkdir build && cd build
cmake -G "Unix Makefiles" ..
cmake --build .
sudo cmake --build . --target install

cd ../../ecos
mkdir build && cd build
cmake ..
make
sudo make install

cd ../../osqp-eigen
mkdir build && cd build
cmake ..
make
sudo make install
