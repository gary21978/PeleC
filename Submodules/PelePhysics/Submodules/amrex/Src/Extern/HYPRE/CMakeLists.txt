add_amrex_define(AMREX_USE_HYPRE NO_LEGACY IF AMReX_HYPRE)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
        continue()
    endif ()

    #
    # This file gets processed if either AMReX_PETSC or AMReX_HYPRE are ON
    # If only AMReX_PETSC is ON, we just need to include this directory
    # and a couple of files to the build.
    #
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    if (NOT AMReX_HYPRE) # only AMReX_PETSC is active

       target_sources(amrex_${D}d
          PRIVATE
          AMReX_Habec_${D}D_K.H
          )

    else ()

       target_sources(amrex_${D}d
          PRIVATE
          AMReX_HypreABecLap.cpp
          AMReX_HypreABecLap2.cpp
          AMReX_HypreABecLap3.cpp
          AMReX_Hypre.cpp
          AMReX_Habec_${D}D_K.H
          AMReX_HypreNodeLap.cpp
          AMReX_HypreIJIface.cpp
          AMReX_HypreMLABecLap.cpp
          AMReX_HypreMLABecLap_${D}D_K.H
          )

    endif ()
endforeach()
