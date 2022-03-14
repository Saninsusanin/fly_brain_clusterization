cmake_minimum_required(VERSION 3.16)
project(metric)

set(CMAKE_CXX_STANDARD 14)

add_executable(metric main.cpp metric.cpp metric.h)

find_package(Boost COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(metric ${Boost_LIBRARIES})
endif()