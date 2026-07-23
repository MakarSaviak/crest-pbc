# This file is part of fmlip-relay.
# SPDX-Identifier: MIT
#
# fmlip-relay: Fortran/C socket bridge for persistent Python ML potential servers.

project(
  'fmlip_relay',
  ['fortran', 'c'],
  version         : '0.1.0',
  license         : 'MIT',
  meson_version   : '>=0.60.0',
  default_options : [
    'buildtype=release',
    'default_library=static',
    'c_std=gnu11',
  ],
)

# ── Compiler flags and static link args ───────────────────────────────────────
# Sets: fc, cc, fc_id, cc_id, static_build
subdir('config')

# Build dir added for generated .mod files
inc_dirs = include_directories('.')

# ═══════════════════════════════════════════════════════════════════════════════
#  OpenMP
# ═══════════════════════════════════════════════════════════════════════════════
# ifx (intel-llvm) triggers a Meson bug in dependency('openmp') because its
# preprocessor output lacks the delimiter strings Meson searches for.
# Work around it by constructing the dependency manually for that compiler.
if fc_id == 'intel-llvm'
  if get_option('openmp')
    omp_dep = declare_dependency(
      compile_args : ['-qopenmp'],
      link_args    : ['-qopenmp'],
    )
  else
    omp_dep = declare_dependency()
  endif
else
  omp_dep = dependency('openmp',
    required : get_option('openmp'),
    language : 'fortran',
  )
endif

_omp_link_dep = omp_dep.found() ? [omp_dep] : []

# ═══════════════════════════════════════════════════════════════════════════════
#  Sources
# ═══════════════════════════════════════════════════════════════════════════════
srcs = []
subdir('src')

# ═══════════════════════════════════════════════════════════════════════════════
#  Library target
# ═══════════════════════════════════════════════════════════════════════════════
fmlip_relay_lib = library(
  meson.project_name(),
  sources             : srcs,
  include_directories : inc_dirs,
  dependencies        : _omp_link_dep,
  install             : not (meson.is_subproject() and get_option('default_library') == 'static'),
)

# ── Exported dependency (what crest links against) ────────────────────────────
fmlip_relay_dep = declare_dependency(
  link_with           : fmlip_relay_lib,
  include_directories : [inc_dirs, fmlip_relay_lib.private_dir_include()],
  dependencies        : _omp_link_dep,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Optional example executable
# ═══════════════════════════════════════════════════════════════════════════════
if get_option('build_example')
  subdir('example')
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Install
# ═══════════════════════════════════════════════════════════════════════════════
install_data(
  'LICENSE',
  install_dir : get_option('datadir') / 'licenses' / meson.project_name(),
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════════════
summary({
  'Fortran compiler' : '@0@ @1@'.format(fc_id, fc.version()),
  'C compiler'       : '@0@ @1@'.format(cc_id, cc.version()),
  'Build type'       : get_option('buildtype'),
  'Static binary'    : static_build,
  'OpenMP'           : omp_dep.found(),
  'Example'          : get_option('build_example'),
}, section : 'fmlip-relay build configuration')
