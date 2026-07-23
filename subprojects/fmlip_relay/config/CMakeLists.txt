# fmlip_relay/config/CMakeLists.txt
# ─────────────────────────────────────────────────────────────────────────────
# Build options, compiler minimum-version checks, and flag configuration.
# Mirrors the gfnff config/ pattern.
# ─────────────────────────────────────────────────────────────────────────────

# ── CMake module path ─────────────────────────────────────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

install(
  DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Fortran module installation subdirectory (compiler-specific to avoid clashes)
set(
  module-dir
  "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)

# ── User options ──────────────────────────────────────────────────────────────
option(WITH_OpenMP       "Enable OpenMP support (parallel server instances)" ON)
option(WITH_TESTS        "Build the test suite"                              OFF)
option(BUILD_EXAMPLE     "Build the example executable"                      OFF)

# ── Default build type ────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type (Debug | Release | RelWithDebInfo | MinSizeRel)" FORCE
  )
  message(STATUS
    "No build type specified – defaulting to '${CMAKE_BUILD_TYPE}'"
  )
endif()

# ── Compiler version requirements ─────────────────────────────────────────────
include("${CMAKE_CURRENT_SOURCE_DIR}/modules/fmlip_relay-utils.cmake")

set(fortran_minimal_versions "GNU;9.0" "Intel;19.0" "IntelLLVM;2021.1")
check_minimal_compiler_version("Fortran" "${fortran_minimal_versions}")

set(c_minimal_versions "GNU;9.0" "Intel;19.0" "IntelLLVM;2021.1")
check_minimal_compiler_version("C" "${c_minimal_versions}")

# ── Compiler flags ────────────────────────────────────────────────────────────
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(fortran_dialect
    "-ffree-line-length-none"
    "-fbacktrace"
  )
  set(fortran_debug_flags
    "-fbounds-check"
    "-ffpe-trap=invalid,zero,overflow"
    "-fcheck=all"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR
       CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  set(fortran_dialect
    "-r8"
    "-align" "array64byte"
    "-traceback"
  )
  set(fortran_debug_flags
    "-check" "all"
    "-fpe0"
  )
else()
  message(WARNING
    "Unrecognised Fortran compiler '${CMAKE_Fortran_COMPILER_ID}' – "
    "no dialect flags set. Compilation may succeed but is untested."
  )
  set(fortran_dialect "")
  set(fortran_debug_flags "")
endif()

# Apply via target compile options rather than global flag strings so they
# don't bleed into C/CXX compilation units.
add_library(fmlip_relay-flags INTERFACE)
target_compile_options(
  fmlip_relay-flags INTERFACE
  "$<$<COMPILE_LANGUAGE:Fortran>:${fortran_dialect}>"
  "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${fortran_debug_flags}>"
)
