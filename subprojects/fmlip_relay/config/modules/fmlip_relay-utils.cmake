# fmlip_relay-utils.cmake
# ─────────────────────────────────────────────────────────────────────────────
# Utility functions for the fmlip_relay CMake build.
# ─────────────────────────────────────────────────────────────────────────────

#
# check_minimal_compiler_version(lang compiler_versions)
#
# Abort with a human-readable error if the active compiler for <lang> is older
# than the minimum required version.
#
# compiler_versions is a flat list of alternating "CompilerID;version" pairs:
#   "GNU;9.0" "Intel;19.0"
#
function(check_minimal_compiler_version lang compiler_versions)
  while(compiler_versions)
    list(POP_FRONT compiler_versions compiler version)
    if("${CMAKE_${lang}_COMPILER_ID}" STREQUAL "${compiler}"
        AND CMAKE_${lang}_COMPILER_VERSION VERSION_LESS "${version}")
      message(FATAL_ERROR
        "${compiler} ${lang} compiler is too old "
        "(found \"${CMAKE_${lang}_COMPILER_VERSION}\", "
        "required >= \"${version}\")"
      )
    endif()
  endwhile()
endfunction()
