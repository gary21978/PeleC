# ---------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the generic SUNDIALS modules
# ---------------------------------------------------------------

# From here we only install the generic SUNDIALS headers. The implementations
# themselves are incorporated in the individual SUNDIALS solver libraries.

install(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add variable sundials_HEADERS with the exported SUNDIALS header files
set(sundials_HEADERS
    sundials_adaptcontroller.h
    sundials_base.hpp
    sundials_context.h
    sundials_context.hpp
    sundials_convertibleto.hpp
    sundials_core.h
    sundials_core.hpp
    sundials_dense.h
    sundials_direct.h
    sundials_errors.h
    sundials_iterative.h
    sundials_linearsolver.h
    sundials_linearsolver.hpp
    sundials_logger.h
    sundials_math.h
    sundials_matrix.h
    sundials_matrix.hpp
    sundials_memory.h
    sundials_memory.hpp
    sundials_nonlinearsolver.h
    sundials_nonlinearsolver.hpp
    sundials_nvector.h
    sundials_nvector.hpp
    sundials_profiler.h
    sundials_profiler.hpp
    sundials_stepper.h
    sundials_types.h
    sundials_version.h)

if(ENABLE_MPI)
  list(APPEND sundials_HEADERS sundials_mpi_errors.h)
endif()

if(ENABLE_CUDA)
  list(APPEND sundials_HEADERS sundials_cuda_policies.hpp)
endif()

if(ENABLE_HIP)
  list(APPEND sundials_HEADERS sundials_hip_policies.hpp)
endif()

if(ENABLE_SYCL)
  list(APPEND sundials_HEADERS sundials_sycl_policies.hpp)
endif()

# Add prefix with complete path to the header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

set(sundials_SOURCES)

if(ENABLE_MPI)
  list(APPEND sundials_SOURCES sundials_mpi_errors.c)
endif()

# Add prefix with complete path to the source files
add_prefix(${SUNDIALS_SOURCE_DIR}/src/sundials/ sundials_SOURCES)

if(ENABLE_MPI)
  set(_link_mpi_if_needed PUBLIC MPI::MPI_C
                          $<$<LINK_LANGUAGE:CXX>:MPI::MPI_CXX>)
endif()

if(SUNDIALS_BUILD_WITH_PROFILING)
  if(ENABLE_CALIPER)
    set(_link_caliper_if_needed PUBLIC caliper)
  endif()
  if(ENABLE_ADIAK)
    set(_link_adiak_if_needed PUBLIC adiak::adiak ${CMAKE_DL_LIBS})
  endif()
endif()

list(LENGTH sundials_SOURCES _sundials_sources_len)
if(_sundials_sources_len EQUAL 0)
  add_library(sundials_core INTERFACE)
  target_link_libraries(sundials_core INTERFACE ${_link_mpi_if_needed})
  target_include_directories(
    sundials_core
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
              $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
              $<BUILD_INTERFACE:${SUNDIALS_SOURCE_DIR}/src/sundials>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  install(TARGETS sundials_core EXPORT sundials-targets)
else()
  # Create a library out of the generic sundials modules
  sundials_add_library(
    sundials_core
    SOURCES ${sundials_SOURCES}
    HEADERS ${sundials_HEADERS}
    INCLUDE_SUBDIR sundials
    LINK_LIBRARIES ${_link_mpi_if_needed}
    OUTPUT_NAME sundials_core
    VERSION ${sundialslib_VERSION}
    SOVERSION ${sundialslib_SOVERSION})
endif()

# Install private headers
install(
  FILES ${SUNDIALS_SOURCE_DIR}/include/sundials/priv/sundials_context_impl.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sundials/priv")

if(NOT _sundials_sources_len EQUAL 0)
  # SUNDIALS_EXPORT macro
  include(GenerateExportHeader)
  generate_export_header(
    sundials_core BASE_NAME SUNDIALS EXPORT_FILE_NAME
    "${PROJECT_BINARY_DIR}/include/sundials/sundials_export.h")
endif()

# F2003 Fortran interface pruned
