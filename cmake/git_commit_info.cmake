# Use git branch and commit hash information in a CMake project
#
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/

get_directory_property( CMDLINEDEFS COMPILE_DEFINITIONS )

if( NOT CMDLINEDEFS MATCHES "${PROJECT_NAME}_GIT_BRANCH" )
	# Get the current working branch
	execute_process(
	  COMMAND git rev-parse --abbrev-ref HEAD
	  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	  OUTPUT_VARIABLE ${PROJECT_NAME}_GIT_BRANCH
	  OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	# make the definition visible to the preprocessor
	add_compile_definitions("${PROJECT_NAME}_GIT_BRANCH=${${PROJECT_NAME}_GIT_BRANCH}")
endif()

if( NOT CMDLINEDEFS MATCHES "${PROJECT_NAME}_GIT_COMMIT_HASH" )
	# Get the latest abbreviated commit hash of the working branch
	execute_process(
	  COMMAND git log -1 --format=%h
	  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	  OUTPUT_VARIABLE ${PROJECT_NAME}_GIT_COMMIT_HASH
	  OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	# make the definition visible to the preprocessor
	add_compile_definitions("${PROJECT_NAME}_GIT_COMMIT_HASH=${${PROJECT_NAME}_GIT_COMMIT_HASH}")
endif()
