add_amrex_define(AMREX_USE_PROFPARSER)
foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
        PRIVATE
        AMReX_BLProfStats.cpp
        AMReX_BLProfUtilities.cpp
        AMReX_CommProfStats.cpp
        AMReX_RegionsProfStats.cpp
        AMReX_AVGDOWN_${D}D.F
        AMReX_BLWritePlotFile.cpp
        AMReX_ProfParserBatch.cpp
    )
endforeach()
