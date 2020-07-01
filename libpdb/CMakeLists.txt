cmake_minimum_required(VERSION 3.10.0)
project(libpdb)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

set(libpdb_SOURCES
	pdb_chars.cpp
	pdb_read.cpp
	pdb_sprntf.cpp
	pdb_sscanf.cpp
	pdb_type.cpp
	pdbinput.cpp
	pdb++.cpp
)

set(libpdb_HEADERS
	pdb++.h
)

add_library(pdb++ ${libpdb_SOURCES} ${libpdb_HEADERS})
# In 2020, all current C++ compilers will surely define the bool type.
# @todo Add a CMake check for whether bool is predefined, and set based on that.
target_compile_definitions(pdb++ PUBLIC -DBOOLPREDEFINED)

# @todo Add configuration information to installed output so those who use
# the library will have things configured correctly based on the target.
