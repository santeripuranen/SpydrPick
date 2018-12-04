# Set default build type if none was specified
set( DEFAULT_BUILD_TYPE "Release" )

if( NOT CMAKE_BUILD_TYPE )
	message( STATUS "Setting build type to \"${DEFAULT_BUILD_TYPE}\" since none was specified by user." )
	set( CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Select build type." FORCE )
endif()
