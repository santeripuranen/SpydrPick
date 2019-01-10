# Version $Id:$

cmake_minimum_required(VERSION 3.1)

#option( ${CMAKE_PROJECT_NAME}_ENABLE_APEGRUNT "Find Apegrunt and, if successful, enable use in ${CMAKE_PROJECT_NAME}" true )

##############
## Apegrunt setup
###

set( ${CMAKE_PROJECT_NAME}_NO_APEGRUNT true CACHE INTERNAL "Don't use Apegrunt, if true" ) # Initialize with default value 

if( ${CMAKE_PROJECT_NAME}_ENABLE_APEGRUNT )
	setup_message( "check for Apegrunt" )
	add_definitions( -DNDEBUG )
	add_subdirectory( apegrunt )

	if( APEGRUNT_FOUND )
		set( ${CMAKE_PROJECT_NAME}_NO_APEGRUNT false CACHE INTERNAL "Don't use Apegrunt, if true" )

		setup_message( "found Apegrunt" )
		setup_message( "include dir: ${APEGRUNT_INCLUDE_DIR}" INDENT )
		setup_message( "library: ${APEGRUNT_LIBRARIES}" INDENT )

		set( APEGRUNT_INCLUDE_DIR ${APEGRUNT_INCLUDE_DIR} CACHE INTERNAL "Apegrunt include directory" )
		set( APEGRUNT_LIBRARIES ${APEGRUNT_LIBRARIES} CACHE INTERNAL "Apegrunt libraries" )
	else()
		setup_message( "WARNING: could not find Apegrunt" )
	endif()		
endif()
if( ${CMAKE_PROJECT_NAME}_NO_APEGRUNT )
	add_definitions( -D${CMAKE_PROJECT_NAME}_NO_APEGRUNT )
	setup_message( "Apegrunt is DISABLED" )
else()
	setup_message( "Apegrunt is enabled" )
endif()	

