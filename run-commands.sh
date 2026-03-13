#!/bin/sh

# Change directory (DO NOT CHANGE!)
repoDir=$(dirname "$(realpath "$0")")
echo $repoDir
cd $repoDir

mkdir -p build
cd build
cmake ..
make -j4
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/tbb_cmake_build/tbb_cmake_build_subdir_release



nsys profile --stats=true ./quicktree -m ../data/matrix_10k.phy -o tree.nwk
# ./quicktree -m ../data/matrix.phy -o tree.nwk

