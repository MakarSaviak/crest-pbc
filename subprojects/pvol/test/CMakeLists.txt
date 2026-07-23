# This file is part of libpvol.
# SPDX-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.17)

###########################################################################
# Fortran (test-drive)
###########################################################################
add_executable(
  "${PROJECT_NAME}-tester"
  coffeine.f90
  lih.f90
  supermol.f90
  timer.f90
  main.f90
  test_pvol.F90
)
target_include_directories(
  "${PROJECT_NAME}-tester"
  PRIVATE "${CMAKE_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "test-drive::test-drive"
)
add_dependencies("${PROJECT_NAME}-tester" "${PROJECT_NAME}-lib")
add_test("${PROJECT_NAME}/pvol" "${PROJECT_NAME}-tester" "pvol")

###########################################################################
# C API
###########################################################################
add_executable("${PROJECT_NAME}-c-tester" main.c)
target_include_directories(
  "${PROJECT_NAME}-c-tester"
  PRIVATE "${CMAKE_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-c-tester"
  PRIVATE "${PROJECT_NAME}-lib"
)
add_test("${PROJECT_NAME}/pvol-C-API" "${PROJECT_NAME}-c-tester")

###########################################################################
# C++ API
###########################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable("${PROJECT_NAME}-cpp-tester" main.cpp)
set_target_properties(
  "${PROJECT_NAME}-cpp-tester"
  PROPERTIES LINKER_LANGUAGE CXX
)
target_include_directories(
  "${PROJECT_NAME}-cpp-tester"
  PRIVATE "${CMAKE_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-cpp-tester"
  PRIVATE "${PROJECT_NAME}-lib"
)
add_test("${PROJECT_NAME}/pvol-C++-API" "${PROJECT_NAME}-cpp-tester")
