# This file is part of libpvol.
# SPDX-Identifier: LGPL-3.0-or-later

# ── Fortran test (test-drive) ──────────────────────────────────────────────────
pvol_tester = executable(
  'pvol-tester',
  sources : files(
    'coffeine.f90',
    'lih.f90',
    'supermol.f90',
    'timer.f90',
    'main.f90',
    'test_pvol.F90',
  ),
  dependencies        : [libpvol_dep, testdrive_dep],
  include_directories : inc_dirs,
  link_language       : 'fortran',
)

test('pvol/pvol', pvol_tester, args : ['pvol'], timeout : 120)

# ── C API test ────────────────────────────────────────────────────────────────
pvol_c_tester = executable(
  'pvol-c-tester',
  sources             : files('main.c'),
  dependencies        : libpvol_dep,
  include_directories : inc_dirs,
)

test('pvol/pvol-C-API', pvol_c_tester, timeout : 120)

# ── C++ API test ──────────────────────────────────────────────────────────────
pvol_cpp_tester = executable(
  'pvol-cpp-tester',
  sources             : files('main.cpp'),
  dependencies        : libpvol_dep,
  include_directories : inc_dirs,
  # No link_language override: g++ drives the link to pull in libstdc++.
  # libgfortran arrives transitively via libpvol_dep.
)

test('pvol/pvol-C++-API', pvol_cpp_tester, timeout : 120)
