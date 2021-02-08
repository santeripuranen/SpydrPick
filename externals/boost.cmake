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
	if( Boost_FOUND ) # is Boost already set up?
		set( ${CMAKE_PROJECT_NAME}_NO_BOOST false CACHE INTERNAL "Don't use Boost, if true" )
	else()
		setup_message( "check for Boost" )
	
		# Use the host's preferred installation prefix, if defined
		if( NOT BOOST_ROOT AND DEFINED ENV{BOOST_ROOT}  )
			set( BOOST_ROOT $ENV{BOOST_ROOT} )
			message( "Got BOOST_ROOT=\"${BOOST_ROOT}\" from shell variable" )
		endif()
		
		if( BOOST_ROOT )
			# turn off system paths if BOOST_ROOT is defined
			set( Boost_NO_SYSTEM_PATHS ON )
		endif()
	
		set( Boost_USE_STATIC_LIBS ON )
		set( Boost_USE_MULTITHREADED TRUE )
		#set( Boost_NO_BOOST_CMAKE ON )
		#set( Boost_DEBUG ON )
	
		find_package( Boost REQUIRED program_options filesystem iostreams timer chrono ) #system #date_time thread
		if( Boost_FOUND )
			set( ${CMAKE_PROJECT_NAME}_NO_BOOST false CACHE INTERNAL "Don't use Boost, if true" )
			setup_message( "found Boost v${Boost_VERSION_STRING}" )
			setup_message( "include dir: ${Boost_INCLUDE_DIRS}" INDENT )
			#setup_message( "library dir: ${Boost_LIBRARY_DIRS}" INDENT )
			setup_message( "libraries: ${Boost_LIBRARIES}" INDENT )
			
			# Work-around for the issue of Boost.CMake targets not being visible in parent/global scope
			foreach( TARGET ${Boost_LIBRARIES} )
				get_target_property(Boost_LOCATION ${TARGET} LOCATION ) # get the full library path
				list( APPEND Boost_LIBRARIES_PATHS ${Boost_LOCATION} )
			endforeach()
			
			# Place into global scope
			set( Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} CACHE INTERNAL "Boost include directory" )
			#set( Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIRS} CACHE INTERNAL "Boost library directory" )
			#set( Boost_LIBRARIES ${Boost_LIBRARIES} CACHE INTERNAL "List of Boost library targets" )
			set( Boost_LIBRARIES ${Boost_LIBRARIES_PATHS} CACHE INTERNAL "List of linkable Boost libraries" )
			#setup_message( "Boost libraries: ${Boost_LIBRARIES}" INDENT )
			
			# stop compiler from nagging about deprecated auto_ptr in boost v1.59.0 and earlier
			set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations" CACHE INTERNAL "" )
		else()
			setup_message( "WARNING: could not find Boost" )
		endif()
	endif()
endif()
if( ${CMAKE_PROJECT_NAME}_NO_BOOST )
	add_definitions( -D${CMAKE_PROJECT_NAME}_NO_BOOST )
	setup_message( "Boost is DISABLED" )
else()
	setup_message( "Boost is enabled" )
endif()
