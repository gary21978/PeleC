foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_AmrCore.cpp
       AMReX_Cluster.cpp
       AMReX_ErrorList.cpp
       AMReX_AmrCore.H
       AMReX_ErrorList.H
       AMReX_FillPatchUtil.cpp
       AMReX_FluxRegister.cpp
       AMReX_FillPatchUtil.H
       AMReX_FillPatchUtil_I.H
       AMReX_FillPatcher.H
       AMReX_InterpBase.H
       AMReX_InterpBase.cpp
       AMReX_MFInterpolater.H
       AMReX_MFInterpolater.cpp
       AMReX_Interpolater.cpp
       AMReX_TagBox.cpp
       AMReX_AmrMesh.cpp
       AMReX_Interpolater.H
       AMReX_TagBox.H
       AMReX_AmrMesh.H
       AMReX_InterpFaceRegister.cpp
       
       )

    if (AMReX_FORTRAN)
       target_sources(amrex_${D}d
          PRIVATE
          AMReX_FillPatchUtil_${D}d.F90
          )
    endif ()

    # These files are needed only if Particles are enabled
    if ( AMReX_PARTICLES )
      target_sources( amrex_${D}d PRIVATE AMReX_AmrParGDB.H AMReX_AmrParticles.H )
    endif ()
endforeach()
