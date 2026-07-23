# fmlip_relay/CMakeLists.txt
# ─────────────────────────────────────────────────────────────────────────────
# Fortran/C socket bridge library for persistent Python ML potential servers.
#
# Quick start:
#   cmake -S . -B build [-DWITH_OpenMP=ON] [-DBUILD_EXAMPLE=ON]
#   cmake --build build
#   cmake --install build --prefix /path/to/install
#
# Consuming in another CMake project:
#   find_package(fmlip_relay REQUIRED)
#   target_link_libraries(myapp PRIVATE fmlip_relay::fmlip_relay)
# ─────────────────────────────────────────────────────────────────────────────

cmake_minimum_required(VERSION 3.17)

project(
  "fmlip_relay"
  VERSION      "0.1.0"
  DESCRIPTION  "Fortran/C bridge for persistent Python ML potential servers"
  LANGUAGES    "Fortran" "C"
)

include(GNUInstallDirs)

# ── Configuration (options, compiler flags, module path) ──────────────────────
set(srcs)
add_subdirectory("config")

# ── OpenMP ────────────────────────────────────────────────────────────────────
# One OpenMP target expected by the library and the example.
# If WITH_OpenMP is OFF we create a dummy interface so all target_link_libraries
# calls below are unconditional and clean.
if(NOT TARGET "OpenMP::OpenMP_Fortran")
  if(WITH_OpenMP)
    find_package("OpenMP" REQUIRED COMPONENTS Fortran)
  else()
    add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
  endif()
endif()

# ── Sources ───────────────────────────────────────────────────────────────────
add_subdirectory("src")

# ── Library target ────────────────────────────────────────────────────────────
add_library("${PROJECT_NAME}-lib" ${srcs})

set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  OUTPUT_NAME                  "${PROJECT_NAME}"
  VERSION                      "${PROJECT_VERSION}"
  SOVERSION                    "${PROJECT_VERSION_MAJOR}"
  POSITION_INDEPENDENT_CODE    TRUE
  WINDOWS_EXPORT_ALL_SYMBOLS   TRUE
  Fortran_MODULE_DIRECTORY     "${PROJECT_BINARY_DIR}/include"
)

target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "OpenMP::OpenMP_Fortran"
  "fmlip_relay-flags"
)

target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

# Ensure the module output directory exists at configure time
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

# ── Install & export ──────────────────────────────────────────────────────────

# Thin interface target that downstream find_package() users link against
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")

install(
  TARGETS
    "${PROJECT_NAME}"
    "${PROJECT_NAME}-lib"
    "fmlip_relay-flags"
  EXPORT  "${PROJECT_NAME}-targets"
  LIBRARY     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME     DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

install(
  EXPORT      "${PROJECT_NAME}-targets"
  NAMESPACE   "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Install generated Fortran .mod files
install(
  DIRECTORY   "${PROJECT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)

# Generate and install CMake package config files so find_package() works
include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fmlip_relay-config.cmake.in"
  "${PROJECT_BINARY_DIR}/cmake/fmlip_relay-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/cmake/fmlip_relay-config-version.cmake"
  VERSION       "${PROJECT_VERSION}"
  COMPATIBILITY "SameMajorVersion"
)

install(
  FILES
    "${PROJECT_BINARY_DIR}/cmake/fmlip_relay-config.cmake"
    "${PROJECT_BINARY_DIR}/cmake/fmlip_relay-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# ── Example executable ────────────────────────────────────────────────────────
if(BUILD_EXAMPLE)
  add_subdirectory("example")
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────
if(WITH_TESTS)
  enable_testing()
  add_subdirectory("tests")
endif()
