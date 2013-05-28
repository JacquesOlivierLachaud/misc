# This is the CMake script for compiling a CGAL application.

project( DGtalCGAL ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()


SET( CMAKE_MODULE_PATH /opt/local/lib/cmake )

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  message(STATUS "Boost_LIBRAIRIES=" ${Boost_LIBRAIRIES})
  message(STATUS "Boost_PROGRAM_OPTIONS_LIBRARY=" ${Boost_PROGRAM_OPTIONS_LIBRARY})
  include_directories( ${Boost_INCLUDE_DIRS} )
  LINK_DIRECTORIES(${Boost_PROGRAM_OPTIONS_LIBRARY})
#  LINK_DIRECTORIES(${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif( Boost_FOUND )



find_package(CGAL COMPONENTS Core )
include( ${CGAL_USE_FILE} )

# find_package(DGtal REQUIRED)
include_directories(${CGAL_INCLUDE_DIR})

FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


SET(SRCs
  simpleExample
  2d-triangulation
  2d-patterns-v2
  3d-patterns
  3d-patterns-v2)

if ( CGAL_FOUND )
  
  FOREACH(FILE ${SRCs})
    add_executable(${FILE} ${FILE})
    target_link_libraries( ${FILE}  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${DGTAL_LIBRARIES} ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  ENDFOREACH(FILE)
  
else()

  message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")

endif()
