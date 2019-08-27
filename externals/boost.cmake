# Version $Id:$

cmake_minimum_required(VERSION 3.1)

#option( ${CMAKE_PROJECT_NAME}_ENABLE_BOOST "Find Boost and, if successful, enable use in ${CMAKE_PROJECT_NAME}" true )

###############
## Boost setup
###

# Make sure the v3 of boost.filesystem library is used
add_definitions( -DBOOST_FILESYSTEM_VERSION=3 )
# When cross-compiling, we need to explicitly state that
# the thread library needs to be linked
add_definitions( -DBOOST_THREAD_USE_LIB )

set( ${CMAKE_PROJECT_NAME}_NO_BOOST true CACHE INTERNAL "Don't use Boost, if true" ) # Initialize with default value
if( ${CMAKE_PROJECT_NAME}_ENABLE_BOOST )
	setup_message( "check for Boost" )
	# List of usable boost versions.

	set( Boost_USE_STATIC_LIBS ON )
	set( Boost_USE_MULTITHREADED TRUE )
	set( Boost_NO_BOOST_CMAKE ON )
	
	find_package( Boost REQUIRED program_options filesystem iostreams timer ) # system # chrono date_time thread
	if( Boost_FOUND )
		set( ${CMAKE_PROJECT_NAME}_NO_BOOST false CACHE INTERNAL "Don't use Boost, if true" )
		setup_message( "found Boost v${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )
		setup_message( "include dir: ${Boost_INCLUDE_DIRS}" INDENT )
		setup_message( "library dir: ${Boost_LIBRARY_DIRS}" INDENT )
		
		# Place into global scope
		set( Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} CACHE INTERNAL "Boost include directory" )
		set( Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIRS} CACHE INTERNAL "Boost library directory" )
		set( Boost_LIBRARIES ${Boost_LIBRARIES} CACHE INTERNAL "List of linkable Boost libraries" )
		
		# stop compiler from nagging about deprecated auto_ptr in boost v1.59.0 and earlier
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations" CACHE INTERNAL "" )
	else()
		setup_message( "WARNING: could not find Boost" )
	endif()
endif()
if( ${CMAKE_PROJECT_NAME}_NO_BOOST )
	add_definitions( -D${CMAKE_PROJECT_NAME}_NO_BOOST )
	setup_message( "Boost is DISABLED" )
else()
	setup_message( "Boost is enabled" )
endif()
