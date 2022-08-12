cmake_minimum_required(VERSION 3.0)
project(unitTest VERSION 0.0.1)

file(GLOB_RECURSE TEST_SRC_LIST *.c)
find_package(OpenMP REQUIRED)
foreach(v ${TEST_SRC_LIST})
  string(REGEX MATCH "test/.*" relative_path ${v})
  string(REGEX REPLACE "test/" "" target_name ${relative_path})
  string(REGEX REPLACE "\\.c" "" target_name ${target_name})

  add_executable(${target_name} ${v})
  add_test(NAME ${target_name} COMMAND ${target_name})
  target_link_libraries(${target_name} sqp)

  if("${target_name}" STREQUAL "test_optimize_sqp_with_blzd")

  endif()

  if(OpenMP_FOUND)
    target_link_libraries(${target_name} OpenMP::OpenMP_C)
  else()

  endif()

endforeach()
