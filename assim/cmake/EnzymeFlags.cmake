# =====================================================================
#  EnzymeFlags.cmake — flags de compilation propres à Enzyme, DÉDIÉS à la
#  couche AD et ACTIVÉS UNIQUEMENT par le mot-clé -DENABLE_ENZYME=ON.
#
#  Aucun de ces flags ne touche le SW2D natif : ils ne s'appliquent qu'à la
#  cible AD (voir CMakeLists.txt) et seulement sous ENABLE_ENZYME.
#
#  Schéma Enzyme (LTO) :
#   - COMPILE : -flto génère du bitcode LLVM (pas d'objet final), sans vectorisation
#     ni déroulage (meilleur pour l'AD). PAS de plugin ici.
#   - LINK    : -flto fusionne le bitcode, lld charge le plugin Enzyme qui génère
#     le code de différenciation sur l'IR fusionné.
# =====================================================================

if(ENABLE_ENZYME)
    set(ENZYME_CLANG_PLUGIN
        "/home/chapumax/Enzyme/build/Enzyme/ClangEnzyme-18.so"
        CACHE STRING "Chemin du plugin Enzyme (clang)")
    set(ENZYME_LLD_PLUGIN
        "/home/chapumax/Enzyme/build/Enzyme/LLDEnzyme-18.so"
        CACHE STRING "Chemin du plugin Enzyme (lld)")

    if(NOT EXISTS ${ENZYME_CLANG_PLUGIN})
        message(FATAL_ERROR "Plugin Enzyme introuvable : ${ENZYME_CLANG_PLUGIN}")
    endif()

    find_program(LLD18 NAMES lld-18 ld.lld-18 lld)
    if(NOT LLD18)
        message(FATAL_ERROR "lld-18 introuvable (sudo apt install lld-18)")
    endif()

    set(ENZYME_COMPILE_FLAGS
        -flto -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops)
    set(ENZYME_LINK_FLAGS
        -flto -fuse-ld=lld-18
        -fpass-plugin=${ENZYME_LLD_PLUGIN}
        -Wl,--load-pass-plugin=${ENZYME_LLD_PLUGIN})

    message(STATUS "Enzyme ACTIVE. plugin=${ENZYME_LLD_PLUGIN}")
    message(STATUS "  ENZYME_COMPILE_FLAGS=${ENZYME_COMPILE_FLAGS}")
    message(STATUS "  ENZYME_LINK_FLAGS=${ENZYME_LINK_FLAGS}")
endif()
