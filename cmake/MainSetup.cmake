
## ###################################################################
## Requirements for C++ stantards
## ###################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## ###################################################################
## Executables destinations
## ###################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin) # general
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CMAKE_BUILD_TYPE} ${PROJECT_BINARY_DIR}/bin) # per-configuration

## ###################################################################
## Get the compiler name
## ###################################################################
get_filename_component(COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME_WE)

## ###################################################################
## Get the version number (for doxygen)
## ###################################################################
file(READ "${CMAKE_SOURCE_DIR}/src/common/version.h" VERSION_CONTENTS)
string(REGEX MATCH "sw2d_version[ \t]*=[ \t]*\"([0-9.]+)\"" _ "${VERSION_CONTENTS}")
set(SW2D_VERSION "${CMAKE_MATCH_1}")
message(STATUS "Detected version for SW2D: ${SW2D_VERSION}")
configure_file(${CMAKE_SOURCE_DIR}/doc/doxyfile ${CMAKE_SOURCE_DIR}/doc/doxyfile_with_version @ONLY)


## ###################################################################
## Linker flags no manifest files 
## ###################################################################
if((COMPILER_NAME STREQUAL "cl") OR (COMPILER_NAME STREQUAL "icx"))
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST:NO")
endif()

## ###################################################################
## Coverage flags
## ###################################################################
if(ENABLE_COVERAGE)
    message(STATUS "Coverage enabled")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g --coverage")
    else()
        message(WARNING "Coverage not supported for compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()
