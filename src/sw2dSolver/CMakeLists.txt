project(sw2dSolver)

## #################################################################
## Sources
## #################################################################
set(${PROJECT_NAME}_HEADERS
  ${COMMON_HEADERS}
  common/swSolverParam.h
  models/common/swBoundaryConditionRatingCurve.h
  models/common/swBoundaryConditionTimeSeries.h
  models/common/swHyp.h
  models/common/swParab.h
  models/common/swSource.h
  models/ddp/swDdp.h
  models/dip/swDip.h
  models/liswes/swLiswes.h
  models/sp/swSp.h
  models/swes/swSwes.h
  models/swParent.h
)

set(${PROJECT_NAME}_SOURCES
  ${COMMON_SOURCES}
  common/swSolverParam.cpp
  main.cpp
  models/common/swBoundaryConditionRatingCurve.cpp
  models/common/swBoundaryConditionTimeSeries.cpp
  models/common/swHyp.cpp
  models/common/swParab.cpp
  models/common/swSource.cpp
  models/common/swSourceFriction.cpp
  models/common/swSourceRunoff.cpp
  models/common/swSourceWind.cpp
  models/swParent.cpp
  models/ddp/swDdp.cpp
  models/ddp/swDdp.cpp
  models/ddp/swDdpHyp.cpp
  models/ddp/swDdpHypGodunov.cpp
  models/ddp/swDdpHypMusclEvr.cpp
  models/ddp/swDdpParab.cpp
  models/ddp/swDdpSourceFriction.cpp
  models/ddp/swDdpSourceRunoffInfiltration.cpp
  models/ddp/swDdpSourceWind.cpp
  models/dip/swDip.cpp
  models/dip/swDipHyp.cpp
  models/dip/swDipHypGodunov.cpp
  models/dip/swDipHypMusclEvr.cpp
  models/dip/swDipParab.cpp
  models/dip/swDipSourceFriction.cpp
  models/dip/swDipSourceRunoffInfiltration.cpp
  models/dip/swDipSourceWind.cpp
  models/liswes/swLiswes.cpp
  models/liswes/swLiswesHyp.cpp
  models/liswes/swLiswesHypGodunov.cpp
  models/liswes/swLiswesHypMusclEvr.cpp
  models/liswes/swLiswesParab.cpp
  models/liswes/swLiswesSourceFriction.cpp
  models/liswes/swLiswesSourceRunoffInfiltration.cpp
  models/liswes/swLiswesSourceWind.cpp
  models/sp/swSp.cpp
  models/sp/swSpHyp.cpp
  models/sp/swSpHypGodunov.cpp
  models/sp/swSpHypMusclEvr.cpp
  models/sp/swSpParab.cpp
  models/sp/swSpSourceFriction.cpp
  models/sp/swSpSourceRunoffInfiltration.cpp
  models/sp/swSpSourceWind.cpp
  models/swes/swSwes.cpp
  models/swes/swSwesHyp.cpp
  models/swes/swSwesHypGodunov.cpp
  models/swes/swSwesHypMusclEvr.cpp
  models/swes/swSwesParab.cpp
  models/swes/swSwesSourceFriction.cpp
  models/swes/swSwesSourceRunoffInfiltration.cpp
  models/swes/swSwesSourceWind.cpp
)

## #################################################################
## Build rules
## #################################################################

add_executable(${PROJECT_NAME}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS} "models/common/swParab.cpp")


#target_link_libraries(${PROJECT_NAME} PRIVATE swCore)

######################################################################
### CMakeLists.txt ends here
