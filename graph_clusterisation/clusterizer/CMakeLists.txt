cmake_minimum_required(VERSION 3.12)
project(hgmc)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_subdirectory(src)
add_subdirectory(logger)
add_subdirectory(app)