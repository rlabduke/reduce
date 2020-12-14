cmake_minimum_required(VERSION 3.10.0)
project(reduce)

# So that we can wrap our library into a shared library using Boost.Python
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(libpdb)
add_subdirectory(toolclasses)
add_subdirectory(reduce_src)