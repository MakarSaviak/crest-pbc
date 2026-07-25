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

printf '%s  %s\n' \
  '02aefedb780a969523dd6c08a8c34999cfef2c4cab33e1d35ec93c52544ec3ce' \
  "$patch_dir/crest-exact-frozen-mtd-gradients.patch" | sha256sum -c -
git apply --check "$patch_dir/crest-exact-frozen-mtd-gradients.patch"
git apply "$patch_dir/crest-exact-frozen-mtd-gradients.patch"

git submodule update --init --recursive subprojects/gfnff
cd "$root/subprojects/gfnff"
printf '%s  %s\n' \
  '81b5a87f1a8fb6ce3e74c625e16d1f193ba90e50b47eca2cb5340775101747dd' \
  "$patch_dir/gfnff-host-guest-fragments.patch" | sha256sum -c -
git apply --check "$patch_dir/gfnff-host-guest-fragments.patch"
git apply "$patch_dir/gfnff-host-guest-fragments.patch"

printf '%s  %s\n' \
  '4bb0cef08cc7de2a382c74d17699992089c252d79c5abe7237824c15905e9756' \
  "$patch_dir/gfnff-exact-frozen-mtd-gradients.patch" | sha256sum -c -
git apply --check "$patch_dir/gfnff-exact-frozen-mtd-gradients.patch"
git apply "$patch_dir/gfnff-exact-frozen-mtd-gradients.patch"

echo "Applied CREST 3.0.2 host-guest workflow, fragment backport, and exact frozen-host MTD acceleration."
