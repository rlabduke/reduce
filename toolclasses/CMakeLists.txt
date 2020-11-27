cmake_minimum_required(VERSION 3.10.0)
project(toolclasses)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

add_definitions(
  -DCHARFUNCMACROS
  -DINCTEMPLATEDEFNS
  -DBOOLPREDEFINED
  -DBRACKETOPERPARMS
  -DLISTFRIENDFIX
)

set(toolclasses_SOURCES
	UseCount.cpp
	Point3d.cpp
	utility.cpp
	column_sprintf.cpp
	column_sscanf.cpp
	DisjointSets.cpp
)

set(toolclasses_HEADERS
	UseCount.h
	Point3d.h
	utility.h
	DisjointSets.h
)

add_library(toolclasses ${toolclasses_SOURCES} ${toolclasses_HEADERS})
# In 2020, all current C++ compilers will surely define the bool type.
# @todo Add a CMake check for whether bool is predefined, and set based on that.
target_compile_definitions(toolclasses PUBLIC
	-DCHARFUNCMACROS -DINCTEMPLATEDEFNS -DBOOLPREDEFINED
	-DBRACKETOPERPARMS -DLISTFRIENDFIX
)
if (NOT WIN32)
  target_link_libraries(toolclasses PUBLIC m)
endif()

# @todo Add configuration information to installed output so those who use
# the library will have things configured correctly based on the target.
