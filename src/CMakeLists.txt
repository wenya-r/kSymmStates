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





option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")

if(ENABLE_UNIT_TESTS)
include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.8.0
)

FetchContent_GetProperties(googletest)

#

if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)

  # Prevent GoogleTest from overriding our compiler/linker options
  # when building with Visual Studio
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  # Prevent GoogleTest from using PThreads
  set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)

  # adds the targers: gtest, gtest_main, gmock, gmock_main
  add_subdirectory(
    ${googletest_SOURCE_DIR}
    ${googletest_BINARY_DIR}
    )

  # Silence std::tr1 warning on MSVC
  if(MSVC)
    foreach(_tgt gtest gtest_main gmock gmock_main)
      target_compile_definitions(${_tgt}
        PRIVATE
          "_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING"
        )
    endforeach()
  endif()
endif()
add_executable(cpp_test "")

target_sources(cpp_test
  PRIVATE
    test.cpp
  )

target_link_libraries(cpp_test
  PRIVATE
    basic
    gtest_main
  )

enable_testing()

add_test(
  NAME google_test
  COMMAND $<TARGET_FILE:cpp_test>
  )








endif()





#add_library(cpp_test test.cpp)
##target_link_libraries(cpp_test basicAlg)
#
#enable_testing()
#
#
##add_test(
##  NAME bash_test
##  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test.sh $<TARGET_FILE:sum_up>
##  )
#
#add_test(
#  NAME cpp_test
#  COMMAND $<TARGET_FILE:cpp_test>
#  )
