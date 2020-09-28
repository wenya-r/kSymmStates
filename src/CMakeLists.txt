cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(overlapK LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(basic 
    STATIC
        basicAlg.hpp
        basicAlg.cpp
)

add_library(basicAlg basicAlg.cpp)
#add_executable(overlap_r createAllStates.cpp)
add_executable(overlap_r rotation.cpp)
target_link_libraries(overlap_r basic)




add_library(cpp_test test.cpp)
#target_link_libraries(cpp_test basicAlg)

enable_testing()


#add_test(
#  NAME bash_test
#  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test.sh $<TARGET_FILE:sum_up>
#  )

add_test(
  NAME cpp_test
  COMMAND $<TARGET_FILE:cpp_test>
  )
