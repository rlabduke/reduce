cmake_minimum_required(VERSION 3.10.0)
project(reduce)

# So that we can wrap our library into a shared library using Boost.Python
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Require C++-11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(libpdb)
add_subdirectory(toolclasses)
add_subdirectory(reduce_src)
