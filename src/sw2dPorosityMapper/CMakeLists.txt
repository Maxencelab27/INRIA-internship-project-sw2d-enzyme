project(sw2dPorosityMapper)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
    ${COMMON_HEADERS}
    swPorosityMapperParam.h
    swPorosityMapperPipeline.h
    swRasterBuildings.h
)

set(${PROJECT_NAME}_SOURCES
    ${COMMON_SOURCES}
    swPorosityMapperParam.cpp
    swPorosityMapperPipeline.cpp
    swRasterBuildings.cpp
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
