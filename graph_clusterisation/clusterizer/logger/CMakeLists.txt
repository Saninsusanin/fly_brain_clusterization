project(logger)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

add_library(logger
        logger.cpp
        logger.h)

find_package(Boost COMPONENTS log log_setup system thread REQUIRED)

if(Boost_FOUND)
    target_link_libraries(logger ${Boost_LIBRARIES} pthread)
endif()