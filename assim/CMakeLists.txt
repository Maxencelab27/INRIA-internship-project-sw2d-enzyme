# =====================================================================
#  CMakeLists.txt — build AUTONOME de la couche d'assimilation / AD (Enzyme).
#
#  - Ne modifie AUCUN fichier natif (référence ../src en entrée).
#  - Les flags Enzyme sont ACTIVÉS UNIQUEMENT par -DENABLE_ENZYME=ON
#    (voir cmake/EnzymeFlags.cmake). Sans ce mot-clé, la cible AD n'est pas
#    construite (elle requiert Enzyme pour le link de __enzyme_autodiff).
#
#  Build :
#    cmake -S assim -B build_assim -DENABLE_ENZYME=ON \
#          -DCMAKE_CXX_COMPILER=clang++-18
#    cmake --build build_assim -j
# =====================================================================
cmake_minimum_required(VERSION 3.16)
project(sw2d_assim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_ENZYME "Build the Enzyme AD assimilation engine" OFF)

set(NATIVE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# Bibliothèques d'optimisation externes (header-only : LBFGSpp, psqn ; .so : optim).
# Par défaut pointées vers le dépôt de développement ; à vendorer dans le release.
set(OPTIM_DIR   "/home/chapumax/sw2d_enzyme/optim"   CACHE PATH "optim (liboptim.so + include)")
set(LBFGSPP_DIR "/home/chapumax/sw2d_enzyme/LBFGSpp" CACHE PATH "LBFGSpp include")
set(PSQN_DIR    "/home/chapumax/sw2d_enzyme/psqn"    CACHE PATH "psqn include")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/EnzymeFlags.cmake)

if(NOT ENABLE_ENZYME)
    message(STATUS "ENABLE_ENZYME=OFF : la cible AD n'est PAS construite.")
    message(STATUS "  -> reconfigurer avec -DENABLE_ENZYME=ON -DCMAKE_CXX_COMPILER=clang++-18")
    return()
endif()

# Sources NATIVES (référencées, NON modifiées). Compilées AVEC la cible AD pour
# qu'Enzyme voie leur bitcode LTO (fonctions natives appelées dans la tape :
# computeMusclSlopes, lambdaMinus/Plus, specificForce...).
file(GLOB_RECURSE NATIVE_COMMON ${NATIVE_ROOT}/common/*.cpp)
file(GLOB_RECURSE NATIVE_SOLVER ${NATIVE_ROOT}/sw2dSolver/*.cpp)
list(FILTER NATIVE_SOLVER EXCLUDE REGEX "/main\\.cpp$")   # pas le main natif

# Sources de la couche AD (noms de fichiers renommés ; noms de FONCTIONS inchangés).
set(AD_SOURCES
    swSpE_core.cpp
    swSpE_hyp.cpp
    swSpE_muscl.cpp
    swSpE_friction.cpp
    swSpE_checkpoint.cpp
    swSpE_ad.cpp
    swSpE_minimize.cpp
    swSpE_gauges.cpp
    swSpE_spAD.cpp
    assimConfig.cpp
    main_enzyme.cpp
)

add_executable(sw2d_enzyme ${AD_SOURCES} ${NATIVE_COMMON} ${NATIVE_SOLVER})

target_include_directories(sw2d_enzyme PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${NATIVE_ROOT}/sw2dSolver/models/sp     # résout les includes relatifs natifs
    ${NATIVE_ROOT}/sw2dSolver/models        #   (../common/..., ../swParent.h) des .cpp AD
    ${OPTIM_DIR}/include
    ${LBFGSPP_DIR}/include
    ${PSQN_DIR}
    ${PSQN_DIR}/inst/include
    /usr/include/eigen3
)
target_compile_definitions(sw2d_enzyme PRIVATE OPTIM_ENABLE_EIGEN_WRAPPERS ENABLE_ENZYME_AD)
target_link_libraries(sw2d_enzyme PRIVATE ${OPTIM_DIR}/liboptim.so)
set_target_properties(sw2d_enzyme PROPERTIES BUILD_RPATH ${OPTIM_DIR})

# Flags Enzyme : appliqués À LA SEULE cible AD, uniquement sous ENABLE_ENZYME.
message(STATUS ">>> Enzyme ACTIVE pour sw2d_enzyme <<<")
target_compile_options(sw2d_enzyme PRIVATE ${ENZYME_COMPILE_FLAGS})
target_link_options(sw2d_enzyme PRIVATE ${ENZYME_LINK_FLAGS})
