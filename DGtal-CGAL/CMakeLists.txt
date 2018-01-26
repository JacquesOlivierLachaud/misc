# This is the CMake script for compiling a CGAL application.

project( DGtalCGAL ) 

cmake_minimum_required (VERSION 2.6) 

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fext-numeric-literals")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(use_cxx11)

use_cxx11()

FIND_PACKAGE(DGtal 0.7 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
# if ( COMMAND cmake_policy )
#   cmake_policy( SET CMP0003 NEW )  
# endif()

# SET( CMAKE_MODULE_PATH /opt/local/lib/cmake ${CMAKE_SOURCE_DIR})

# FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
# if ( Boost_FOUND )
#   message(STATUS "Boost and boost_program_options found.")
#   message(STATUS "Boost_LIBRAIRIES=" ${Boost_LIBRAIRIES})
#   message(STATUS "Boost_PROGRAM_OPTIONS_LIBRARY=" ${Boost_PROGRAM_OPTIONS_LIBRARY})
#   include_directories( ${Boost_INCLUDE_DIRS} )
#   LINK_DIRECTORIES(${Boost_PROGRAM_OPTIONS_LIBRARY})
# #  LINK_DIRECTORIES(${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
# endif( Boost_FOUND )


 # -------------------------------------------------------------------------
 # This test is for instance used for ITK v3.x. ITK forces a limited
 # template depth. We remove this option.
 # --------------------------------------------------------------------------
 if (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")
    message( "Warning: some package has enabled a limited template depth for the C++ compiler." )
    message( "         Disabling option -ftemplate-depth-xx in CMAKE_CXX_FLAGS." )
    set( CMAKE_CXX_FLAGS_TMP ${CMAKE_CXX_FLAGS} )
    STRING( REGEX REPLACE "-ftemplate-depth-[0-9]*" ""
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_TMP}" )
     message ("         CMAKE_CXX_FLAGS=" ${CMAKE_CXX_FLAGS} )
  endif (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")

# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(DGtalTools_VERSION_MAJOR 0)
SET(DGtalTools_VERSION_MINOR 8)
SET(DGtalTools_VERSION_PATCH beta)
SET(DGTALTOOLS_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(PROJECT_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(VERSION ${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH})

SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )




# # DGtal
# FIND_PACKAGE(DGtal REQUIRED)
# INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
# LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# # Cairo
# FIND_PACKAGE(Cairo REQUIRED)
# IF(CAIRO_FOUND)
#     INCLUDE_DIRECTORIES(${CAIRO_INCLUDE_DIRS})
#     LINK_DIRECTORIES(${CAIRO_LIBRAIRIES})
#     message(STATUS "Cairo found in INCLUDE=" ${CAIRO_INCLUDE_DIRS})
#     message(STATUS "Cairo found in LIBRARY=" ${CAIRO_LIBRAIRIES})
# ELSE(CAIRO_FOUND)
#     message(FATAL_ERROR "cairo not found. Check the cmake variables associated to this package or disable it." )
# ENDIF(CAIRO_FOUND)


# SET(SRCs
# #  simpleExample
#   2d-triangulation
#   2d-triangulation-f
#   2d-triangulation-color
#   3d-triangulation
#   3d-triangulation-2
#   3d-triangulation-3
#   2d-patterns-v2
#   3d-patterns
#   3d-patterns-v2
#   inside-delaunay-3d )

# if ( CGAL_FOUND )
  
#   FOREACH(FILE ${SRCs})
# if (  WITH_CGAL )
#   SET( CGAL_TESTS_SRC
#         3d-patterns-v2
#         3d-triangulation
#         plane-analysis )
#   FOREACH(FILE ${CGAL_TESTS_SRC})
#     add_executable(${FILE} ${FILE})
#     target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
#   ENDFOREACH(FILE)
# else ( WITH_CGAL )
#   message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")
# endif ( WITH_CGAL )

find_package(CGAL COMPONENTS Core )
include( ${CGAL_USE_FILE} )
include_directories(${CGAL_INCLUDE_DIR})

# DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# Cairo

FIND_PATH(CAIRO_INCLUDE_DIRS cairo.h
    $ENV{CAIRODIR}/include
    $ENV{CAIRODIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/include
    /usr/include
    /usr/include/cairo # MT
    /usr/local/include/cairo #brew, manual
    /sw/include # Fink
    /opt/local/include # DarwinPorts
    /opt/local/include/cairo # DarwinPorts # MT
    /opt/csw/include # Blastwave
    /opt/include/cairo
    /usr/X11R6/include/cairo
)


find_library(CAIRO_LIBRAIRIES
             NAMES cairo
             PATHS /usr/lib
                  /usr/local/lib
                 ENV CAIRO_ROOT
                 ENV LD_LIBRARY_PATH
                 ENV LIBRARY_PATH
             PATH_SUFFIXES cairo
            )

SET(CAIRO_FOUND 0)
IF(CAIRO_LIBRAIRIES AND CAIRO_INCLUDE_DIRS)
  SET(CAIRO_FOUND 1)
ENDIF(CAIRO_LIBRAIRIES AND CAIRO_INCLUDE_DIRS)

#FIND_PACKAGE(Cairo REQUIRED)
IF(CAIRO_FOUND)
    INCLUDE_DIRECTORIES(${CAIRO_INCLUDE_DIRS})
    LINK_DIRECTORIES(${CAIRO_LIBRAIRIES})
    message(STATUS "Cairo found in INCLUDE=" ${CAIRO_INCLUDE_DIRS})
    message(STATUS "Cairo found in LIBRARY=" ${CAIRO_LIBRAIRIES})
ELSE(CAIRO_FOUND)
    message(FATAL_ERROR "cairo not found. Check the cmake variables associated to this package or disable it." )
ENDIF(CAIRO_FOUND)


SET(SRCs
  simpleExample
  2d-triangulation
  2d-triangulation-f
  2d-triangulation-color
  tv-triangulation-color
  tv-image
  3d-triangulation
  3d-triangulation-2
  3d-triangulation-3
  2d-patterns-v2
  3d-patterns
  3d-patterns-v2)

if ( CGAL_FOUND )
  
  FOREACH(FILE ${SRCs})
    add_executable(${FILE} ${FILE})
    target_link_libraries( ${FILE}  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CAIRO_LIBRAIRIES} ${DGTAL_LIBRARIES} ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  ENDFOREACH(FILE)
  
else()

  message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")

endif()
