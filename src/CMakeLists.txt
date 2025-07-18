# Gpuspline sources
# https://github.com/gpufit/Gpuspline

if( MSVC ) # link runtime statically with MSVC
  foreach( type ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE} )
    string( TOUPPER ${type} TYPE )
    foreach( flags CMAKE_C_FLAGS_${TYPE} CMAKE_CXX_FLAGS_${TYPE} )
      get_property( help CACHE ${flags} PROPERTY HELPSTRING )
      string( REPLACE "/MD" "/MT" ${flags} "${${flags}}" )
      set( ${flags} "${${flags}}" CACHE STRING "${help}" FORCE )
    endforeach()
  endforeach()
endif()

function( add_launcher target executable arguments working_directory )
  if( MSVC12 OR MSVC14 )
    file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/${target}.vcxproj.user
"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
"<Project ToolsVersion=\"14.0\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">\n"
"  <PropertyGroup>\n"
"    <LocalDebuggerCommand>${executable}</LocalDebuggerCommand>\n"
"    <LocalDebuggerCommandArguments>${arguments}</LocalDebuggerCommandArguments>\n"
"    <LocalDebuggerWorkingDirectory>${working_directory}</LocalDebuggerWorkingDirectory>\n"
"  </PropertyGroup>\n"
"</Project>\n"
    )
  endif()
endfunction()

# Default values

# Set the default value for the precision
set( DEFAULT_DOUBLE_PRECISION OFF )

# Set the default value for the Matlab working directory
set( DEFAULT_Matlab_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/src/matlab/examples" )

# Precision

set( DOUBLE_PRECISION ${DEFAULT_DOUBLE_PRECISION} CACHE BOOL "Build with double precision")

if( DEFINED DOUBLE_PRECISION AND DOUBLE_PRECISION)
    add_definitions( -DSPLINE_DOUBLE )
endif()

# MATLAB

find_package( Matlab COMPONENTS MX_LIBRARY )
if( Matlab_FOUND )
  find_program( Matlab_EXECUTABLE matlab
    PATHS "${Matlab_ROOT_DIR}/bin" PATH_SUFFIXES win32 win64 NO_DEFAULT_PATH )
  function( add_matlab_launcher target )
    set( paths "${CMAKE_BINARY_DIR}/$(Configuration)" "${Matlab_WORKING_DIRECTORY}" )
    list( GET paths -1 working_directory )
    string( REPLACE ";" "','" paths "${paths}" )
    set( arguments "-r addpath('${paths}');" )
    add_launcher( ${target} "${Matlab_EXECUTABLE}" "${arguments}" "${working_directory}" )
  endfunction()
endif()

# MATLAB Launcher

if( Matlab_FOUND )
    file(REMOVE ${CMAKE_BINARY_DIR}/Gpuspline.sln)
    if ( "${Matlab_WORKING_DIRECTORY}" STREQUAL "" )
        set( Matlab_WORKING_DIRECTORY ${DEFAULT_Matlab_WORKING_DIRECTORY} CACHE PATH "MATLAB working directory" FORCE)
    endif()
    set( Matlab_WORKING_DIRECTORY ${DEFAULT_Matlab_WORKING_DIRECTORY} CACHE PATH "MATLAB working directory" )
    add_custom_target( RUN_MATLAB )
    set_property( TARGET RUN_MATLAB PROPERTY FOLDER CMakePredefinedTargets )
    add_dependencies( RUN_MATLAB MATLAB_SPLINE_PACKAGE )
    add_matlab_launcher( RUN_MATLAB )
endif()

# Python

find_package( PythonInterp )
if( PYTHONINTERP_FOUND )
  function( add_python_launcher target )
    set( paths "${CMAKE_BINARY_DIR}/$(Configuration)" ${ARGN} )
    list( GET paths -1 working_directory )
    string( REPLACE ";" "')\nsys.path.append('" paths "${paths}" )
    set( arguments "-i -c \"import sys\nsys.path.append('${paths}')\"" )
    add_launcher( ${target} "${PYTHON_EXECUTABLE}" "${arguments}" "${working_directory}" )
  endfunction()
endif()

# Documentation (Latex->PDF)

find_package( LATEX COMPONENTS PDFLATEX )
if( LATEX_PDFLATEX_FOUND )
  set( latex_manual_directory "${CMAKE_SOURCE_DIR}/docs/_build/latex" )
  # run latex two times to get the references right
  add_custom_target( DOCUMENTATION_PDFLATEX
    COMMAND ${CMAKE_COMMAND} -E
      chdir ${latex_manual_directory} "${PDFLATEX_COMPILER}" -interaction=nonstopmode Gpuspline.tex
    COMMAND ${CMAKE_COMMAND} -E
      chdir ${latex_manual_directory} "${PDFLATEX_COMPILER}" -interaction=nonstopmode Gpuspline.tex
    COMMENT "Converting documentation to PDF"
  )
  set_property( TARGET DOCUMENTATION_PDFLATEX PROPERTY FOLDER CMakePredefinedTargets )
else()
  message( WARNING "PDFLatex NOT found - cannot convert documentation to PDF!" )
endif()

# splines

set(splinesHeaders
    spline.h
    spline_classes.h
    definitions.h
    equation_system.h
    math_utils.h
    bspline_1d.h
	bspline_fast_cubic_basis_evaluate.h
	natural_bspline_1d.h
	natural_bspline_nd.h
	lsq_bspline_1d.h
	libs/fitpack/gpuspline_fitpack_functions.h
)

set( splinesSources
    spline.def
    spline_1d.cpp
    spline_2d.cpp
    spline_3d.cpp
	spline_4d.cpp
	spline_5d.cpp
    spline_1d_c_interface.cpp
    spline_2d_c_interface.cpp
    spline_3d_c_interface.cpp
	spline_4d_c_interface.cpp
	spline_5d_c_interface.cpp
	bspline_1d.cpp
	bspline_5d.cpp
	natural_bspline_1d.cpp
	natural_bspline_nd.cpp
	natural_bspline_nd_c_interface.cpp
    lsq_bspline_1d.cpp
    equation_system.cpp
	regrid_data.cpp
	libs/fitpack/gpuspline_fitpack_functions.cpp
)

add_library( splines SHARED
    ${splinesHeaders}
    ${splinesSources}
)

target_include_directories(splines
    PUBLIC ${PROJECT_SOURCE_DIR}
)

target_compile_definitions(splines PRIVATE SPLINES_BUILD_DLL)

set_target_properties(splines PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    CXX_VISIBILITY_PRESET hidden
)

# === spline_tests executable ===

add_executable( splines_tests
    spline_tests.cpp
)

target_include_directories( splines_tests PRIVATE ${PROJECT_SOURCE_DIR} )
target_link_libraries( splines_tests splines )

set_target_properties( splines_tests
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        CXX_VISIBILITY_PRESET hidden
)

# === bspline_tests executable ===

add_executable(bspline_tests
    bspline_tests.cpp
)

target_link_libraries(bspline_tests
    PRIVATE splines
)

target_include_directories(bspline_tests
    PRIVATE ${PROJECT_SOURCE_DIR}
)

set_target_properties(bspline_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    CXX_VISIBILITY_PRESET hidden
)


# === natural_bspline_tests executable ===

add_executable(natural_bspline_tests
    natural_bspline_tests.cpp
)

target_link_libraries(natural_bspline_tests
    PRIVATE splines
)

target_include_directories(natural_bspline_tests
    PRIVATE ${PROJECT_SOURCE_DIR}
)

set_target_properties(natural_bspline_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    CXX_VISIBILITY_PRESET hidden
)



# === lsq_bspline_tests executable ===

add_executable(lsq_bspline_tests
    lsq_bspline_tests.cpp
)

target_link_libraries(lsq_bspline_tests
    PRIVATE splines
)

target_include_directories(lsq_bspline_tests
    PRIVATE ${PROJECT_SOURCE_DIR}
)

set_target_properties(lsq_bspline_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    CXX_VISIBILITY_PRESET hidden
)

# Matlab binding

add_subdirectory( matlab )

# Python binding

add_subdirectory( python )

