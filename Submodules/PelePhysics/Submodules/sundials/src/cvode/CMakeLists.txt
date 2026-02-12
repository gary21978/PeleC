# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban, and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODE\n\")")

# Add variable cvode_SOURCES with the sources for the CVODE library
set(cvode_SOURCES)

# Add variable cvode_HEADERS with the exported CVODE header files
set(cvode_HEADERS cvode.h   
                  cvode_ls.h cvode_proj.h)

# Add prefix with complete path to the CVODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/cvode/ cvode_HEADERS)

# Build fused kernel libraries
if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(cvode_fused_gpu.cpp PROPERTIES LANGUAGE CUDA)

    sundials_add_library(
      sundials_cvode_fused_cuda
      SOURCES cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS PRIVATE USE_CUDA
      LINK_LIBRARIES PUBLIC sundials_core PRIVATE sundials_nveccuda
      OUTPUT_NAME sundials_cvode_fused_cuda
      VERSION ${cvodelib_VERSION}
      SOVERSION ${cvodelib_SOVERSION})
  endif()

  if(BUILD_NVECTOR_HIP)
    sundials_add_library(
      sundials_cvode_fused_hip
      SOURCES cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS PRIVATE USE_HIP
      LINK_LIBRARIES PUBLIC sundials_core PRIVATE sundials_nvechip
      OUTPUT_NAME sundials_cvode_fused_hip
      VERSION ${cvodelib_VERSION}
      SOVERSION ${cvodelib_SOVERSION})
  endif()

  set(_fused_link_lib)
endif()

list(LENGTH cvode_SOURCES _cvode_sources_len)
if(_cvode_sources_len EQUAL 0)
  add_library(sundials_cvode INTERFACE)
  target_link_libraries(sundials_cvode INTERFACE sundials_core ${_fused_link_lib})
  target_include_directories(
    sundials_cvode
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
              $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
              $<BUILD_INTERFACE:${SUNDIALS_SOURCE_DIR}/src/sundials>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  install(TARGETS sundials_cvode EXPORT sundials-targets)
else()
  # Create the library
  sundials_add_library(
    sundials_cvode
    SOURCES ${cvode_SOURCES}
    HEADERS ${cvode_HEADERS}
    INCLUDE_SUBDIR cvode
    LINK_LIBRARIES PUBLIC sundials_core
    LINK_LIBRARIES # Link to stubs so examples work.
                   PRIVATE ${_fused_link_lib}
    OUTPUT_NAME sundials_cvode
    VERSION ${cvodelib_VERSION}
    SOVERSION ${cvodelib_SOVERSION})
endif()

# Finished CVODE
message(STATUS "Added CVODE module")

# F2003 Fortran interface pruned
