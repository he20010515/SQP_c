cmake_minimum_required(VERSION 3.16)
project(sqp_c)
include_directories ("include")

aux_source_directory("src/function" DIR_SRCS)
aux_source_directory("src/linarg" DIR_SRCS)
aux_source_directory("src/util" DIR_SRCS)
aux_source_directory("src/optimize" DIR_SRCS)
LINK_LIBRARIES(m)

add_library(sqp STATIC ${DIR_SRCS})
# add_executable(sqp ${DIR_SRCS})
enable_testing()
add_subdirectory(test)