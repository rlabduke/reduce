# This is designed to be build as a subdirectory build where the libpdb and toolclasses
# targets have already been defined.

cmake_minimum_required(VERSION 3.10.0)
project(reduce_src)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

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
	RotMethyl.cpp
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
	RotMethyl.h
	RotAromMethyl.h
	RotDonor.h
  FlipMemo.h
	CliqueList.h
	AtomDescr.h
	PDBrec.h
	MoveableNode.h
	hybrid_36_c.h
	neighbors.h
)

if (WIN32)
	set(DICT_HOME "C:/usr/local/reduce_wwPDB_het_dict.txt")
	set(DICT_OLD "C:/usr/local/reduce_het_dict.txt")
else ()
	set(DICT_HOME "/usr/local/reduce_wwPDB_het_dict.txt")
	set(DICT_OLD "/usr/local/reduce_het_dict.txt")
endif ()

add_executable(reduce ${reduce_SOURCES} ${reduce_HEADERS})
# In 2020, all current C++ compilers will surely define the bool type.
# @todo Add a CMake check for whether bool is predefined, and set based on that.
target_compile_definitions(reduce PUBLIC
	-DAROMATICS_ACCEPT_HBONDS -DCHARFUNCMACROS -DBOOLPREDEFINED
	-DBRACKETOPERPARMS -DLEFT_JUSTIFY_NUC_RES_OK
  -DLISTFRIENDFIX -DINCTEMPLATEDEFNS
	-DHET_DICTIONARY="${DICT_HOME}"
	-DHET_DICTOLD="${DICT_OLD}"
)
target_include_directories(reduce PRIVATE ../libpdb ../toolclasses)
target_link_libraries(reduce pdb++ toolclasses)

# @todo Add configuration information to installed output so those who use
# the library will have things configured correctly based on the target.
