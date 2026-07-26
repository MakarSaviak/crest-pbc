#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch="$root/patches/gfnff-exact-eeq-workspace-20260725.patch"

cd "$root"
printf '%s  %s\n' \
  'fe9850da1aa99e202f1383c8d607beb63824ef0a4285ae5b79d8ce251e4bd018' \
  "$patch" | sha256sum -c -

git apply --check "$patch"
git apply "$patch"

echo "Applied persistent exact EEQ work arrays on top of the exact frozen-host cache checkpoint."
