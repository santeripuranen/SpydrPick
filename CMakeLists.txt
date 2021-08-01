# Version $Id: $

cmake_minimum_required(VERSION 3.3)

# Disable in-source builds to prevent source tree corruption
include( cmake/disable_insource_builds.cmake NO_POLICY_SCOPE )

# For some reason CMake (at least 3.1.1, 3.3.1) does not always honor the
# compiler environment variables unless we define them in the local scope
set( CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER} )
set( CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER} )

project( SPYDRPICK )
set( CMAKE_CXX_STANDARD 17 )
set( ${PROJECT_NAME}_INCLUDE_DIR ${${PROJECT_NAME}_SOURCE_DIR}/include )

include( cmake/SpydrPick_options.cmake NO_POLICY_SCOPE )
include( cmake/git_commit_info.cmake NO_POLICY_SCOPE )
include( cmake/set_default_build_type.cmake NO_POLICY_SCOPE )

add_subdirectory( externals )

add_subdirectory( src )
#add_subdirectory( test )
