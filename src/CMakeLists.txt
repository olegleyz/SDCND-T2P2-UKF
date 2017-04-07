cmake_minimum_required(VERSION 3.6)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp tracking.cpp tracking.h ukf.cpp ukf.h measurement_package.h )
add_executable(src ${SOURCE_FILES})