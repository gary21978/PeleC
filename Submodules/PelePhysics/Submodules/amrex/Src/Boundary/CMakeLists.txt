foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_Mask.cpp
       AMReX_MultiMask.cpp
       AMReX_LO_BCTYPES.cpp
       AMReX_YAFluxRegister.H
       AMReX_YAFluxRegister_K.H
       AMReX_YAFluxRegister_${D}D_K.H
       )

    if (NOT (D EQUAL 1))
       target_sources(amrex_${D}d
          PRIVATE
          AMReX_EdgeFluxRegister.cpp
          )
    endif ()

    if (AMReX_FORTRAN)
       target_sources(amrex_${D}d
          PRIVATE
          )
    endif ()
endforeach()
