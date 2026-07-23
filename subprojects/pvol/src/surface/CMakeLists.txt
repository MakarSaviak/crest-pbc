# This file is part of libpvol.
# SPDX-Identifier: LGPL-3.0-or-later
#
# NOTE: sasa.f90, surface.f90, tesspoints.f90 exist in this directory but are
# intentionally not included here — they were absent from the original build
# systems as well.  Add them when ready.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/vdwrad.f90"
  "${dir}/lebedev.f90"
  "${dir}/bisectsearch.f90"
  "${dir}/math_wrapper.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
