# Version $Id: $

# spydrpick src

############################
## Add sources and includes
###

include_directories(
	${SPYDRPICK_INCLUDE_DIR}
	${APEGRUNT_INCLUDE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
	${Boost_INCLUDE_DIR}
	${TBB_INCLUDE_DIRS}
)

link_directories( ${Boost_LIBRARY_DIRS} )

set( SPYDRPICK_SOURCES
	SpydrPick.cpp
	SpydrPick_options.cpp
) # *.cpp *.hpp *.cc

set( ARACNE_SOURCES
	aracne.cpp
	aracne_options.cpp
)

#################################
## Add libraries and executables
###

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# general optimization flags	
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -O3 -mavx -fabi-version=6 -static-libstdc++" )
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -ftree-vectorize -fwhole-program" )
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -flto -ffat-lto-objects" ) # -flto-report" )

# release build flags
set( SPYDRPICK_GCC_RELEASE_FLAGS "-w -DNDEBUG" )
#set( SPYDRPICK_GCC_RELEASE_FLAGS "${SPYDRPICK_GCC_RELEASE_FLAGS} -fdata-sections -ffunction-sections -Wl,--gc-sections" ) # surprisingly, results in bigger binary, but no perceivable effect on performance
set( SPYDRPICK_GCC_RELEASE_FLAGS "${SPYDRPICK_GCC_RELEASE_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden -Wl,--strip-all" )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} ${SPYDRPICK_GCC_RELEASE_FLAGS}")

# debug build flags
set( SPYDRPICK_GCC_DEBUG_FLAGS "-pg -g -ftree-vectorizer-verbose=2" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} ${SPYDRPICK_GCC_DEBUG_FLAGS}")

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wignored-attributes" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS}" )

set( CMAKE_LINKER "ld.gold" )

if( SPYDRPICK_SOURCES )
	add_executable( SpydrPick
		${SPYDRPICK_SOURCES}
	)
	target_link_libraries( SpydrPick apegrunt )
	set_target_properties( SpydrPick PROPERTIES COMPILE_FLAGS "--std=c++14 ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS}" )
endif()

if( ARACNE_SOURCES )
	add_executable( aracne
		${ARACNE_SOURCES}
	)
	set_target_properties( aracne PROPERTIES COMPILE_FLAGS "--std=c++14 ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS}" )
endif()

# Add Boost libraries
if( NOT SPYDRPICK_NO_BOOST )
	target_link_libraries( SpydrPick ${Boost_LIBRARIES} )
	target_link_libraries( aracne ${Boost_LIBRARIES} )
endif()

# Add TBB libraries
if( NOT SPYDRPICK_NO_TBB )
	target_link_libraries( SpydrPick ${TBB_LIBRARIES} )
	target_link_libraries( aracne ${TBB_LIBRARIES} )
endif()

# Add pthreads library if building in UNIX
if( UNIX )
	target_link_libraries( SpydrPick pthread )
	target_link_libraries( aracne pthread )
endif( UNIX )

# Prevent linking against shared libraries on OS X;
# Apple gcc always links against a shared version of a library if present,
# regardless of -Bstatic or equivalent linker flags.
if(APPLE)
	set_target_properties( SpydrPick PROPERTIES LINK_SEARCH_END_STATIC TRUE )
	set_target_properties( aracne PROPERTIES LINK_SEARCH_END_STATIC TRUE )
endif(APPLE)
