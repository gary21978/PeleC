#
# -DAMREX_USE_SENSEI_INSITU
#   will compile adaptors and bridge and enables the instrumentation in
#   amrex::Amr
#
# -DAMREX_NO_SENSEI_AMR_INST
#   will disable instrumentation in amrex::Amr in that case explicit
#   instrumentation must be available in the simulation code
#
add_amrex_define( AMREX_USE_SENSEI_INSITU IF AMReX_SENSEI NO_LEGACY )
add_amrex_define( AMREX_NO_SENSEI_AMR_INST IF AMReX_NO_SENSEI_AMR_INST NO_LEGACY )

#
# Mesh based adaptors and bridges
#
foreach(D IN LISTS AMReX_SPACEDIM)
    set ( amrex_sensei_sources
       AMReX_AmrMeshDataAdaptor.cpp
       AMReX_AmrMeshInSituBridge.cpp
       AMReX_InSituBridge.cpp
       AMReX_InSituUtils.cpp )

    if( AMReX_AMRLEVEL )
       list ( APPEND amrex_sensei_sources
          AMReX_AmrDataAdaptor.cpp
          AMReX_AmrInSituBridge.cpp )
    endif()

    #
    # Particle based adaptors and bridges
    #
    if ( AMReX_PARTICLES )
       list ( APPEND amrex_sensei_sources

       if( AMReX_AMRLEVEL )
          list ( APPEND amrex_sensei_sources
       endif()
    endif ()

    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )
    target_sources(amrex_${D}d PRIVATE ${amrex_sensei_sources} )
endforeach()
