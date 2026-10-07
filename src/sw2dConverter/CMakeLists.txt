project(sw2dConverter)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
    ${COMMON_HEADERS}
    swConverterParam.h
)

set(${PROJECT_NAME}_SOURCES
    ${COMMON_SOURCES}
    swConverterParam.cpp
    main.cpp
)

## #################################################################
## Build rules
## #################################################################

add_executable(${PROJECT_NAME}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS})


#target_link_libraries(${PROJECT_NAME} PRIVATE swCore)

######################################################################
### CMakeLists.txt ends here
