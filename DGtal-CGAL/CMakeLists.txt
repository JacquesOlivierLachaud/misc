# This is the CMake script for compiling a CGAL application.

project( DGtalCGAL ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()


find_package(CGAL COMPONENTS Core )
include( ${CGAL_USE_FILE} )

# find_package(DGtal REQUIRED)
include_directories(${CGAL_INCLUDE_DIR})

FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


SET(SRCs
  simpleExample)

if ( CGAL_FOUND )
  
  FOREACH(FILE ${SRCs})
    add_executable(${FILE} ${FILE})
    target_link_libraries( ${FILE}  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${DGTAL_LIBRARIES} )
  ENDFOREACH(FILE)
  
else()

  message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")

endif()
