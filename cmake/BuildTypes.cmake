## ###################################################################
## Default build type (Release) and compilation options (Warning All, Optimize 3)
## ###################################################################

# Set default build type to Release if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Get the compiler name
get_filename_component(COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME_WE)

# Set compilation flags for all build types
if((COMPILER_NAME STREQUAL "cl") OR (COMPILER_NAME STREQUAL "icx"))
  # "cl" = MSVC and "icx" = Intel (both for Windows)
  set(CMAKE_CXX_FLAGS "/EHsc")
  set(CMAKE_CXX_FLAGS_DEBUG "/MTd /Od /Zi /RTC1 /Wall")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "/MT /O1 /DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "/MT /Ox /DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MT /Ox /Zi /DNDEBUG")
elseif((COMPILER_NAME STREQUAL "c++") OR (COMPILER_NAME STREQUAL "g++"))
  # Both c++ and g++ correspond to the gcc C++ compiler (under Windows and MacOs/Linux)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "")
  else()
    set(CMAKE_CXX_FLAGS "-static")
  endif()
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
else()
  message(WARNING "Unknown compiler specified : ${CMAKE_CXX_COMPILER}. CMake flags will be set on default values.") 
endif()

message(STATUS CMAKE_BUILD_TYPE: "${CMAKE_BUILD_TYPE}")
message(STATUS CMAKE_CXX_FLAGS: "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")