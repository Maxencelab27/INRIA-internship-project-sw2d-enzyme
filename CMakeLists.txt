## ###################################################################
# With this project() call, variables $PROJECT_BINARY_DIR and $PROJECT_SOURCE_DIR are automatically created
# $PROJECT_BINARY_DIR corresponds to the directory from which cmake is run
# $PROJECT_SOURCE_DIR corresponds to the directory provided as an argument to the cmake command
## ###################################################################
cmake_minimum_required(VERSION 3.10.0)
project(sw2d) # Main project name

## ###################################################################
## cmake modules folder # (where cmakes searches infos)
## ###################################################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

## ###################################################################
## Define build types
## ###################################################################
include(BuildTypes)

## ###################################################################
## SW2D main setup
## ###################################################################
include(MainSetup) # Some setups that should not be modified

## ###################################################################
## Input (where there are C++ files to compile)
## ###################################################################
add_subdirectory(src)
add_subdirectory(test/unit)

### CMakeLists.txt ends here