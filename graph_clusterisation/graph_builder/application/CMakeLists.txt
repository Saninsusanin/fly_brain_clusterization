cmake_minimum_required(VERSION 3.16)
project(application)

set(CMAKE_CXX_STANDARD 14)

add_executable(application main.cpp graph_builder.cpp graph_builder.h distances.cpp kernels.cpp kernels.h distances.h)


find_package(Boost COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(application ${Boost_LIBRARIES})
endif()