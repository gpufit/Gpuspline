# Gpuspline
# https://github.com/gpufit/Gpuspline

# CMake

cmake_minimum_required( VERSION 3.11 )
set_property( GLOBAL PROPERTY USE_FOLDERS ON )

if( NOT PROJECT_NAME )
  project( Gpuspline VERSION 1.0.0 )
  include( CTest )
endif()

if( NOT CMAKE_CXX_STANDARD )
	set( CMAKE_CXX_STANDARD 14 )
endif()

add_subdirectory( src )