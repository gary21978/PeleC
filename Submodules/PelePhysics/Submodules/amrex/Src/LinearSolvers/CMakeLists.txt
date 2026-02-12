foreach(D IN LISTS AMReX_SPACEDIM)
    #
    # Sources in subdirectory MLMG
    #
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/MLMG>)

    target_sources(amrex_${D}d
       PRIVATE
       MLMG/AMReX_MLMG_${D}D_K.H
       MLMG/AMReX_MLNodeLinOp_${D}D_K.H
       MLMG/AMReX_MLNodeLinOp.cpp
       MLMG/AMReX_MLCellABecLap_${D}D_K.H
       MLMG/AMReX_MLABecLap_${D}D_K.H
       MLMG/AMReX_MLALap_${D}D_K.H
       MLMG/AMReX_MLPoisson_${D}D_K.H
       )

    if (D EQUAL 3)
       target_sources(amrex_${D}d
          PRIVATE
          )
    endif ()

    if (AMReX_LINEAR_SOLVERS_EM)
       target_sources(amrex_${D}d
          PRIVATE
          )

       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLEBNodeFDLaplacian.cpp
          MLMG/AMReX_MLEBNodeFDLap_${D}D_K.H
          MLMG/AMReX_MLNodeTensorLaplacian.cpp
          MLMG/AMReX_MLNodeTensorLap_${D}D_K.H
       )
    endif ()

    if (AMReX_LINEAR_SOLVERS_INCFLO)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLNodeABecLap_${D}D_K.H
          MLMG/AMReX_MLNodeLaplacian.cpp
          MLMG/AMReX_MLNodeLaplacian_sync.cpp
          MLMG/AMReX_MLNodeLaplacian_sten.cpp
          MLMG/AMReX_MLNodeLaplacian_misc.cpp
          MLMG/AMReX_MLNodeLap_${D}D_K.H
          MLMG/AMReX_MLTensorOp.cpp
          MLMG/AMReX_MLTensor_${D}D_K.H
       )
    endif ()

    if (AMReX_EB AND NOT D EQUAL 1)
       target_sources(amrex_${D}d
          PRIVATE
          MLMG/AMReX_MLEBABecLap.cpp
          MLMG/AMReX_MLEBABecLap_F.cpp
          MLMG/AMReX_MLEBABecLap_${D}D_K.H
          )

       if (AMReX_LINEAR_SOLVERS_INCFLO)
          target_sources(amrex_${D}d
             PRIVATE
             MLMG/AMReX_MLEBTensorOp.cpp
             MLMG/AMReX_MLEBTensor_${D}D_K.H
             )
       endif ()
    endif ()

    if (AMReX_FORTRAN)
       target_sources(amrex_${D}d
          PRIVATE
          )
    endif ()

    if (AMReX_HYPRE AND AMReX_LINEAR_SOLVERS_INCFLO)
       target_sources(amrex_${D}d
          PRIVATE
          )
    endif ()

    if (D EQUAL 3)

       target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/OpenBC>)

       target_sources(amrex_${D}d
          PRIVATE
          OpenBC/AMReX_OpenBC.cpp
          )
    endif ()
endforeach()
