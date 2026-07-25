#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch_dir="$root/patches"
combined="${TMPDIR:-/tmp}/gfnff-exact-frozen-caches-20260725.patch"

cd "$root"
cat "$patch_dir"/gfnff-exact-frozen-caches-20260725.patch.part-* > "$combined"
printf '%s  %s\n' \
  'db97fc4ee03f682509a72c2ded3646e804882639db1ab872587d36050e506a27' \
  "$combined" | sha256sum -c -

git apply --check "$combined"
git apply "$combined"

echo "Applied exact frozen-host GFN-FF caches on top of the active-force-loop checkpoint."
