# This file is part of libpvol.
# SPDX-Identifier: LGPL-3.0-or-later
#
# libpvol is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libpvol is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with libpvol.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.17)

project(
  "pvol"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "0.0.1"
  DESCRIPTION "Hydrostatic pressure force-field library"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")

###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################

# OpenMP
if(NOT TARGET "OpenMP::OpenMP_Fortran")
  if(WITH_OpenMP)
    find_package("OpenMP" REQUIRED)
    list(APPEND libs OpenMP::OpenMP_Fortran)
  else()
    add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
  endif()
endif()

# LAPACK / BLAS
if(NOT TARGET "LAPACK::LAPACK")
  find_package(LAPACK REQUIRED)
endif()
list(
  APPEND libs
  LAPACK::LAPACK
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
)

# test-drive
if(NOT TARGET "test-drive::test-drive" AND WITH_TESTS)
  find_package("test-drive" REQUIRED)
endif()

# pvol sources
set(srcs)
add_subdirectory("src")

###############################################################################
############################### Static Library ################################
###############################################################################
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${libs}"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
##################### Install and Export Instructions #########################
###############################################################################
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT "${PROJECT_NAME}-targets"
  NAMESPACE "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
install(
  FILES "COPYING" "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

###############################################################################
############################## Standalone app #################################
###############################################################################
if(build_exe)
  add_subdirectory("app")
endif()

###############################################################################
############################### Python bindings ###############################
###############################################################################
if(PYTHON_BINDINGS)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 REQUIRED)
  python_add_library(_libpvol src/pylibpvol.cpp WITH_SOABI)
  target_link_libraries(_libpvol
    PRIVATE
    "${PROJECT_NAME}-lib"
    pybind11::headers
  )
  target_compile_definitions(_libpvol PRIVATE VERSION_INFO=${PROJECT_VERSION})
  install(TARGETS _libpvol DESTINATION pylibpvol)
endif()

###############################################################################
################################## Testing ####################################
###############################################################################
if(WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()
