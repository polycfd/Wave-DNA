#!/bin/tcsh
rm -r CMakeCache.txt CMakeFiles Makefile WaveDNA cmake_install.cmake
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Debug
make
