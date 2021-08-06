# This is designed to be build as a subdirectory build where the libpdb and toolclasses
# targets have already been defined.

cmake_minimum_required(VERSION 3.10.0)
project(reduce_src)

set(USE_PYTHON TRUE)
FIND_PACKAGE(PythonInterp)
if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python)
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  else()
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  endif()
else()
    set(USE_PYTHON FALSE)
    message("Python not found")
endif()
if (NOT Boost_FOUND)
    set(USE_PYTHON FALSE)
endif()

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else ()
  add_definitions(
    -DHET_DICTIONARY="/usr/local/reduce_wwPDB_het_dict.txt"
    -DHET_DICTOLD="/usr/local/reduce_het_dict.txt"
  )
endif ()

add_definitions(
  -DAROMATICS_ACCEPT_HBONDS
  -DCHARFUNCMACROS
  -DBOOLPREDEFINED
  -DBRACKETOPERPARMS
  -DLEFT_JUSTIFY_NUC_RES_OK
  -DLISTFRIENDFIX
  -DINCTEMPLATEDEFNS
)

set(reduce_SOURCES
  GraphToHoldScores.cpp
  reduce.cpp
  CTab.cpp
  ElementInfo.cpp
  StdResH.cpp
  ResBlk.cpp
  AtomConn.cpp
  AtomPositions.cpp
  DotSph.cpp
  Mover.cpp
  Rot.cpp
  Rot3Fold.cpp
  RotAromMethyl.cpp
  RotDonor.cpp
  FlipMemo.cpp
  CliqueList.cpp
  AtomDescr.cpp
  PDBrec.cpp
  MoveableNode.cpp
  hybrid_36_c.c
)

set(reduce_HEADERS
  GraphToHoldScores.h
  CTab.h
  ElementInfo.h
  StdResH.h
  ResBlk.h
  AtomConn.h
  AtomPositions.h
  DotSph.h
  Mover.h
  Rot.h
  Rot3Fold.h
  RotAromMethyl.h
  RotDonor.h
  FlipMemo.h
  CliqueList.h
  AtomDescr.h
  PDBrec.h
  MoveableNode.h
  hybrid_36_c.h
  reduce.h
	neighbors.h
)

if (WIN32)
	set(DICT_HOME "C:/usr/local/reduce_wwPDB_het_dict.txt")
	set(DICT_OLD "C:/usr/local/reduce_het_dict.txt")
else ()
	set(DICT_HOME "/usr/local/reduce_wwPDB_het_dict.txt")
	set(DICT_OLD "/usr/local/reduce_het_dict.txt")
endif ()

add_library(reducelib ${reduce_SOURCES} ${reduce_HEADERS})
# In 2020, all current C++ compilers will surely define the bool type.
# @todo Add a CMake check for whether bool is predefined, and set based on that.
target_compile_definitions(reducelib PUBLIC
	-DAROMATICS_ACCEPT_HBONDS -DCHARFUNCMACROS -DBOOLPREDEFINED
	-DBRACKETOPERPARMS -DLEFT_JUSTIFY_NUC_RES_OK
  -DLISTFRIENDFIX -DINCTEMPLATEDEFNS
	-DHET_DICTIONARY="${DICT_HOME}"
	-DHET_DICTOLD="${DICT_OLD}"
)
target_include_directories(reducelib PUBLIC ../libpdb ../toolclasses)
target_link_libraries(reducelib PUBLIC pdb++ toolclasses)

if (USE_PYTHON)
  message(STATUS "Building Boost.Python library")
  message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
  message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
  message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
  message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
  include_directories(../libpdb ../toolclasses ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  # Note: To make this work on Windows, you need to copy the boost_python*.dll file into the
  # executable directory to avoid an error about "DLL load failed: The specified module could not be found"
  PYTHON_ADD_MODULE(mmtbx_reduce_ext reduce_bpl.cpp)
  target_link_libraries(mmtbx_reduce_ext PRIVATE reducelib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

  add_custom_command(
    TARGET mmtbx_reduce_ext POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/reduce_src/reduce.py
      ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

add_executable(reduce main.cpp)
target_include_directories(reduce PUBLIC ./)
target_link_libraries(reduce reducelib)

install(TARGETS reduce DESTINATION bin)
install(TARGETS reduce DESTINATION bin)

# @todo Add configuration information to installed output so those who use
# the library will have things configured correctly based on the target.
