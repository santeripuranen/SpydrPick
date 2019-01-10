# Version $Id:$

cmake_minimum_required(VERSION 3.1)

#option( ${CMAKE_PROJECT_NAME}_ENABLE_TBB "Find TBB and, if successful, enable use in ${CMAKE_PROJECT_NAME}" true )

##############
## Threading Building Blocks setup
###

set( ${CMAKE_PROJECT_NAME}_NO_TBB true CACHE INTERNAL "Don't use TBB, if true" ) # Initialize with default value 
if( ${CMAKE_PROJECT_NAME}_ENABLE_TBB )
	setup_message( "check for TBB" )
	# If FindTBB.cmake is not present in your system, then
	# get it from https://github.com/Kitware/VTK/blob/master/CMake/FindTBB.cmake
	# FindTBB.cmake can be installed anywhere as long as
	# the CMAKE_MODULE_PATH (shell) environment variable
	# is set to point to that location.
	find_package( TBB QUIET )
	if( TBB_FOUND )
		set( ${CMAKE_PROJECT_NAME}_NO_TBB false CACHE INTERNAL "Don't use TBB, if true" )
		setup_message( "found TBB interface v${TBB_INTERFACE_VERSION}" )
		setup_message( "include dirs: ${TBB_INCLUDE_DIRS}" INDENT )
		setup_message( "libraries: ${TBB_LIBRARIES}" INDENT )

		set( TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIRS} CACHE INTERNAL "TBB include directories" )
		set( TBB_LIBRARIES ${TBB_LIBRARIES} CACHE INTERNAL "TBB libraries" )
	else()
		setup_message( "WARNING: could not find TBB headers and libraries" )
	endif()
endif()
if( ${CMAKE_PROJECT_NAME}_NO_TBB )
	add_definitions( -D${CMAKE_PROJECT_NAME}_NO_TBB )
	setup_message( "TBB is DISABLED" )
else()
	setup_message( "TBB is enabled" )
endif()
