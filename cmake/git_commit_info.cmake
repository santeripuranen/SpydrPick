# Use git branch and commit hash information in a CMake project
#
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions("-D${PROJECT_NAME}_GIT_COMMIT_HASH=${${PROJECT_NAME}_GIT_COMMIT_HASH}")
add_definitions("-D${PROJECT_NAME}_GIT_BRANCH=${${PROJECT_NAME}_GIT_BRANCH}")
