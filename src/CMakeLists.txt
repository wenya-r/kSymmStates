cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(overlapK LANGUAGES CXX)


add_library(basic 
    STATIC
        basicAlg.hpp
        basicAlg.cpp
)


#add_executable(overlap_r createAllStates.cpp)
add_executable(overlap_r rotation.cpp)
target_link_libraries(overlap_r basic)

