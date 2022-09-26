#!/bin/bash
sudo apt install libconfig++-dev -y
cd /usr/include
sudo ln -s eigen3/Eigen Eigen

if [ ! -d "~/third_party" ] ; then
    echo "Directory ~/third_party DOES NOT exists."
    mkdir ~/third_party
fi

cd ~/third_party
echo "Installing OSQP..."
git clone --recursive https://github.com/osqp/osqp
cd osqp
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .
sudo cmake --build . --target install


cd ~/third_party
echo "Installing ECOS..."
git clone --recursive https://github.com/vkotaru/ecos.git
cd ecos
mkdir build
cd build
cmake ..
make
sudo make install

cd ~/third_party
echo "Installing OSQP::Eigen..."
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build
cd build
cmake ..
make
sudo make install

sudo ldconfig
