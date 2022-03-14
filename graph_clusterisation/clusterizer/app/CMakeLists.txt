project(app)

include_directories(../src)
include_directories(../logger)
add_executable(run_app main.cpp command_line.cpp command_line.h)

find_package(Boost COMPONENTS random program_options log log_setup system filesystem thread REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(run_app src logger ${Boost_LIBRARIES} pthread)
endif()