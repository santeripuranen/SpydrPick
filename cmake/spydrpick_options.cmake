# Version $Id:$

cmake_minimum_required(VERSION 3.3)

option( ${PROJECT_NAME}_ENABLE_APEGRUNT "Find Apegrunt and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_BOOST "Find Boost and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_TBB "Find TBB and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_COMPILER_INTRINSICS "Find compiler intrinsics headers and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_GPROF "Generate instrumented binaries for profiling with gprof" false ) # off by default
option( ${PROJECT_NAME}_ENABLE_DOXYGEN "Find Doxygen and enable documentation generation" false ) 
