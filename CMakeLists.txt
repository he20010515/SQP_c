cmake_minimum_required(VERSION 3.16)
project(sqp_c)
include_directories ("include")

aux_source_directory("src/function" DIR_SRCS)
aux_source_directory("src/linarg" DIR_SRCS)
aux_source_directory("src/util" DIR_SRCS)
aux_source_directory("src/log" DIR_SRCS)
aux_source_directory("src/optimize" DIR_SRCS)

add_library(sqp STATIC ${DIR_SRCS})

LINK_LIBRARIES(m)

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "found openmp") 
    # set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS}, ${OPENMP_C_FLAGS})
    # # set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}, ${OPENMP_CXX_FLAGS})
    target_link_libraries(sqp OpenMP::OpenMP_C) 
else()
    message(FATAL_ERROR "openmp not found!")
endif()


# add_executable(sqp ${DIR_SRCS})
enable_testing()
add_subdirectory(test)