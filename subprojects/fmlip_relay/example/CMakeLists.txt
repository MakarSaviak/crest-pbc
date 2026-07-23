# fmlip_relay/example/CMakeLists.txt
# ─────────────────────────────────────────────────────────────────────────────
# Builds the example executable.
#
# When built as part of the main project (BUILD_EXAMPLE=ON), it links
# directly against the in-tree fmlip_relay-lib target.
#
# To build the example standalone against an installed fmlip_relay:
#
#   cmake -S . -B build \
#         -DCMAKE_PREFIX_PATH=/path/to/fmlip_relay/install \
#         -DSTANDALONE_EXAMPLE=ON
#   cmake --build build
# ─────────────────────────────────────────────────────────────────────────────

if(DEFINED STANDALONE_EXAMPLE)
  # ── Standalone build ───────────────────────────────────────────────────────
  cmake_minimum_required(VERSION 3.17)
  project("fmlip_relay-example" LANGUAGES "Fortran" "C")
  find_package("fmlip_relay" REQUIRED)
  find_package("OpenMP")   # optional – don't fail if absent

  add_executable("run_example" "example_usage.f90")
  target_link_libraries(
    "run_example"
    PRIVATE
    "fmlip_relay::fmlip_relay"
    "$<$<TARGET_EXISTS:OpenMP::OpenMP_Fortran>:OpenMP::OpenMP_Fortran>"
  )

else()
  # ── In-tree build ──────────────────────────────────────────────────────────
  add_executable("run_example" "example_usage.f90")
  target_link_libraries(
    "run_example"
    PRIVATE
    "${PROJECT_NAME}-lib"
    "$<$<TARGET_EXISTS:OpenMP::OpenMP_Fortran>:OpenMP::OpenMP_Fortran>"
  )
  # Make sure Fortran .mod files are visible
  target_include_directories(
    "run_example"
    PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/../include"
    "${PROJECT_BINARY_DIR}/include"
  )
endif()
