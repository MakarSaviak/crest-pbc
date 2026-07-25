#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch_dir="$root/patches/crest-3.0.2-host-guest"
combined="${TMPDIR:-/tmp}/crest-3.0.2-host-guest.patch"

cd "$root"
cat "$patch_dir"/crest-3.0.2-host-guest.patch.part-* > "$combined"
printf '%s  %s\n' \
  'cdfc64821556e633cea61da9ff2e2a1265b395c9ee7ddffe52d49a94f14e94eb' \
  "$combined" | sha256sum -c -

git apply --check "$combined"
git apply "$combined"

git submodule update --init --recursive subprojects/gfnff
cd "$root/subprojects/gfnff"
printf '%s  %s\n' \
  '81b5a87f1a8fb6ce3e74c625e16d1f193ba90e50b47eca2cb5340775101747dd' \
  "$patch_dir/gfnff-host-guest-fragments.patch" | sha256sum -c -
git apply --check "$patch_dir/gfnff-host-guest-fragments.patch"
git apply "$patch_dir/gfnff-host-guest-fragments.patch"

echo "Applied CREST 3.0.2 host–guest workflow and GFN-FF fragment backport."
