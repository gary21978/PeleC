# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# src level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Add all of the shared SUNDIALS components
add_subdirectory(sundials)
add_subdirectory(nvector)
add_subdirectory(sunmatrix)
add_subdirectory(sunlinsol)
add_subdirectory(sunmemory)

# ARKODE library
if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()
