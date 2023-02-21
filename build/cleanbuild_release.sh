#!/bin/tcsh
rm -r CMakeCache.txt CMakeFiles Makefile WaveDNA cmake_install.cmake cmake.run.macosx.x86_64
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
make
