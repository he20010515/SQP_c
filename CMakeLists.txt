cmake_minimum_required(VERSION 3.20)
project(sqp_c)
include_directories ("include")
aux_source_directory("src" DIR_SRCS)

aux_source_directory("src/function" DIR_SRCS)
aux_source_directory("src/linarg" DIR_SRCS)
aux_source_directory("src/util" DIR_SRCS)


add_executable(sqp ${DIR_SRCS})