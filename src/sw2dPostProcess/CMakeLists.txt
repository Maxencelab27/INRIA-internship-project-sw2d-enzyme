project(sw2dPostProcess)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
    ${COMMON_HEADERS}
    swPostProcessParam.h
    swPostProcessPipeline.h
)

set(${PROJECT_NAME}_SOURCES
    ${COMMON_SOURCES}
    swPostProcessParam.cpp
    swPostProcessPipeline.cpp
    main.cpp
)

## #################################################################
## Build rules
## #################################################################

add_executable(${PROJECT_NAME}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS})

######################################################################
### CMakeLists.txt ends here
