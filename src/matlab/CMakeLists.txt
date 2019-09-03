
# Matlab spline binding

if( NOT Matlab_FOUND )
  message( STATUS "No suitable Matlab and/or MX_Library version found - skipping spline Matlab binding!" )
  return()
endif()

# Matlab MEX FILES

add_library( spline_coefficientsMex SHARED
  spline_coefficientsMex.cpp
  )
  
add_library( spline_interpolateMex SHARED
  spline_interpolateMex.cpp
  )
  
  add_library( spline_valuesMex SHARED
  spline_valuesMex.cpp
  )
  
set_target_properties( spline_coefficientsMex
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

set_target_properties( spline_interpolateMex
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

set_target_properties( spline_valuesMex
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

target_include_directories( spline_coefficientsMex PRIVATE ${Matlab_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
target_link_libraries( spline_coefficientsMex splines ${Matlab_LIBRARIES} )

target_include_directories( spline_interpolateMex PRIVATE ${Matlab_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
target_link_libraries( spline_interpolateMex splines ${Matlab_LIBRARIES} )

target_include_directories( spline_valuesMex PRIVATE ${Matlab_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
target_link_libraries( spline_valuesMex splines ${Matlab_LIBRARIES} )

if( WIN32 )
  SET(CMAKE_SHARED_LINKER_FLAGS "/export:mexFunction")
endif()

add_matlab_launcher( spline_coefficientsMex "${CMAKE_CURRENT_SOURCE_DIR}" )
add_matlab_launcher( spline_interpolateMex "${CMAKE_CURRENT_SOURCE_DIR}" )
add_matlab_launcher( spline_valuesMex "${CMAKE_CURRENT_SOURCE_DIR}" )

# Matlab splines package
set( build_directory "${CMAKE_BINARY_DIR}/${CMAKE_CFG_INTDIR}/matlab" )
set( package_files
  "${CMAKE_CURRENT_SOURCE_DIR}/spline_coefficients.m"
  "${CMAKE_CURRENT_SOURCE_DIR}/spline_values.m"
  "${CMAKE_CURRENT_SOURCE_DIR}/spline_interpolate.m"
)

set( binary_spline $<TARGET_FILE:splines> )
set( binary_spline_mex $<TARGET_FILE:spline_coefficientsMex> )
set( binary_interpolation_mex $<TARGET_FILE:spline_interpolateMex> )
set( binary_spline_values_mex $<TARGET_FILE:spline_valuesMex> )

add_custom_target( MATLAB_SPLINE_PACKAGE ALL
    COMMAND ${CMAKE_COMMAND} -E
        make_directory ${build_directory}
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${package_files} ${build_directory}
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${binary_spline} ${build_directory}
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${binary_spline_mex} ${build_directory}
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${binary_interpolation_mex} ${build_directory}
    COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${binary_spline_values_mex} ${build_directory}
    COMMENT "Creating splines Matlab package"
)

set_property( TARGET MATLAB_SPLINE_PACKAGE PROPERTY FOLDER CMakePredefinedTargets )
add_dependencies( MATLAB_SPLINE_PACKAGE  splines spline_coefficientsMex spline_interpolateMex spline_valuesMex)