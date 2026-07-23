# fmlip_relay/src/CMakeLists.txt
# ─────────────────────────────────────────────────────────────────────────────
# Collects all library sources and appends them to the `srcs` list in the
# parent scope (same pattern as gfnff).
# ─────────────────────────────────────────────────────────────────────────────

list(
  APPEND srcs
  "${CMAKE_CURRENT_SOURCE_DIR}/c/socket_utils.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/mlip_client.f90"
)
set(srcs "${srcs}" PARENT_SCOPE)
