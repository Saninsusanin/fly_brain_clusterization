project(src)

add_library(src
        algorithms/genetic.cpp
        algorithms/leiden.cpp
        utils/utils.h
        algorithms/genetic_algorithm.h
        algorithms/leiden_algorithm.h
        utils/utils.cpp
        algorithms/louvain.cpp
        algorithms/louvain_algorithm.h
        data_structures/nodes.cpp
        data_structures/nodes.h
        data_structures/partition.cpp
        algorithms/common.cpp
        algorithms/common.h
        data_structures/graph.cpp
        data_structures/graph.h
        data_structures/community.cpp
        data_structures/community.h
        data_structures/partition.h
        data_structures/genetic_config.h)
include_directories(../logger)

find_package(Boost COMPONENTS random log log_setup system thread REQUIRED)

if(Boost_FOUND)
    target_link_libraries(src logger ${Boost_LIBRARIES} pthread)
endif()