cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES tpsalib.h mathfunc.h mathfunc.cpp polymap.h polymap.cpp rational.h, nlmap.h)

add_library(tpsalib SOURCE_FILES)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

