#!/usr/bin/env bash
set -euo pipefail
umask 077
export LC_ALL=C

if (( $# != 1 )); then
  printf 'usage: %s SOURCE_ROOT\n' "$0" >&2
  exit 64
fi

source_root=$(readlink -f -- "$1")
test -d "$source_root"
test ! -e "$source_root/.git"

write_file_manifest() (
  cd "$source_root"
  find . -type f -not -path './.git' -not -path './.git/*' -print0 \
    | sort -z | xargs -0 -r sha256sum
)

write_symlink_manifest() (
  cd "$source_root"
  find . -type l -not -path './.git' -not -path './.git/*' \
    -printf '%p -> %l\n' | sort
)

{
  printf 'CREST_SOURCE_TREE_MANIFEST_V1_FILES\n'
  write_file_manifest
  printf 'CREST_SOURCE_TREE_MANIFEST_V1_SYMLINKS\n'
  write_symlink_manifest
} | sha256sum | awk '{print $1}'
