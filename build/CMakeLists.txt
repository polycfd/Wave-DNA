cmake_minimum_required(VERSION 3.12)

project (WaveDNA)
set(CMAKE_C_COMPILER /usr/bin/gcc)

include_directories($ENV{WaveDNA_DIR}/src/include)
FILE (GLOB_RECURSE MYFILES ABSOLUTE  $ENV{WaveDNA_DIR}/src/*.c)

set (mylibs m)

foreach(arg ${myincludes})
  IF (arg MATCHES "-I")
    STRING(REGEX REPLACE "-I" "" myinc ${arg})
    message("Additional include: ${myinc}")
    include_directories(${myinc})
  ENDIF(arg MATCHES "-I")
endforeach(arg ${myincludes})

foreach(arg ${MYLIBS})
  STRING(REGEX REPLACE "lib" "" myl1 ${arg})
	STRING(REGEX REPLACE ".a$" "" myl2 ${myl1})
	set(mylibs ${myl2} ${mylibs} ${myl2})
endforeach(arg)

if(NOT CMAKE_BUILD_TYPE)
   message(STATUS "Optimisation: No optimisation specified.")
   set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
   message(STATUS "Optimization: Debug")
   set(CMAKE_C_FLAGS_DEBUG "-Wall -g")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
   message(STATUS "Optimization: Release")
   add_definitions(-DNDEBUG)
   set(CMAKE_C_FLAGS_RELEASE "-Wall -Werror -O3")
endif()


set(mylibdirs ${mylibdirs} /usr/lib64/)
link_directories(${mylibdirs})
add_executable(WaveDNA ${MYFILES})
target_link_libraries(WaveDNA ${mylibs})
