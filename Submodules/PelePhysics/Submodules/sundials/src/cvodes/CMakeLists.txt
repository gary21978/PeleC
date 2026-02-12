# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the CVODES library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODES\n\")")

# Add variable cvodes_SOURCES with the sources for the CVODES library
set(cvodes_SOURCES
    cvodea.c
    cvodea_io.c
    cvodes_io.c
    cvodes_nls.c
    cvodes_nls_sim.c
    cvodes_nls_stg.c
  cvodes_nls_stg1.c)

# Add variable cvodes_HEADERS with the exported CVODES header files
set(cvodes_HEADERS)

# Add prefix with complete path to the CVODES header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# Create the library
sundials_add_library(
  sundials_cvodes
  SOURCES ${cvodes_SOURCES}
  HEADERS ${cvodes_HEADERS}
  INCLUDE_SUBDIR cvodes
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
    sundials_sunmemsys_obj
    sundials_nvecserial_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
    sundials_sunnonlinsolnewton_obj
    sundials_sunnonlinsolfixedpoint_obj
  OUTPUT_NAME sundials_cvodes
  VERSION ${cvodeslib_VERSION}
  SOVERSION ${cvodeslib_SOVERSION})

# Finished CVODES
message(STATUS "Added CVODES module")

# F2003 Fortran interface pruned
