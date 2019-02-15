# Version $Id: $

# SpydrPick src

############################
## Add sources and includes
###

include_directories(
	${SPYDRPICK_INCLUDE_DIR}
	${APEGRUNT_INCLUDE_DIR}
	$(ARACNE_INCLUDE_DIR)
	${CMAKE_CURRENT_BINARY_DIR}
	${Boost_INCLUDE_DIR}
	${TBB_INCLUDE_DIRS}
)

link_directories( ${Boost_LIBRARY_DIRS} )

set( SPYDRPICK_SOURCES
	SpydrPick.cpp
	SpydrPick_options.cpp
	ARACNE_options.cpp
) # *.cpp *.hpp *.cc

#################################
## Add libraries and executables
###

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# general optimization flags	
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -O3 -mavx2 -static-libstdc++" )
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -ftree-vectorize -fwhole-program" )
set( SPYDRPICK_GCC_OPTIMIZATION_FLAGS "${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} -flto -ffat-lto-objects" ) # -flto-report" )

# release build flags
set( SPYDRPICK_GCC_RELEASE_FLAGS "-w -DNDEBUG -fvisibility=hidden -fvisibility-inlines-hidden -Wl,--strip-all" )
set( SPYDRPICK_CLANG_RELEASE_FLAGS "-w -DNDEBUG -fvisibility=hidden -fvisibility-inlines-hidden -Wl,-s" )

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} ${SPYDRPICK_GCC_RELEASE_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} ${SPYDRPICK_CLANG_RELEASE_FLAGS}")
endif()

# debug build flags
set( SPYDRPICK_GCC_DEBUG_FLAGS "-pg -g -ftree-vectorizer-verbose=2" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS} ${SPYDRPICK_GCC_DEBUG_FLAGS}")

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wignored-attributes" )
#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SPYDRPICK_GCC_OPTIMIZATION_FLAGS}" )

# set preferred linker
set( CMAKE_LINKER "ld.gold" )

if( SPYDRPICK_SOURCES )
	add_executable( SpydrPick
		${SPYDRPICK_SOURCES}
	)
	target_link_libraries( SpydrPick libapegrunt )
	#target_link_libraries( SpydrPick ${APEGRUNT_LIBRARIES} )
	set_target_properties( SpydrPick PROPERTIES COMPILE_FLAGS "--std=c++14" )
endif()

# Add Boost libraries
if( NOT SPYDRPICK_NO_BOOST )
	target_link_libraries( SpydrPick ${Boost_LIBRARIES} )
endif()

# Add TBB libraries
if( NOT SPYDRPICK_NO_TBB )
	target_link_libraries( SpydrPick ${TBB_LIBRARIES} )
endif()

# Add rt and pthreads library if building in UNIX
if( UNIX )
	target_link_libraries( SpydrPick pthread rt )
endif( UNIX )

# Prevent linking against shared libraries on OS X;
# Apple gcc always links against a shared version of a library if present,
# regardless of -Bstatic or equivalent linker flags.
if(APPLE)
	set_target_properties( SpydrPick PROPERTIES LINK_SEARCH_END_STATIC TRUE )
endif(APPLE)
