#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch="$root/patches/gfnff-exact-eeq-host-block-20260725.patch"

cd "$root"
printf '%s  %s\n' \
  'b3f196ec284761add3e874ac7630f37f99c1cb0cd068f4c789a8e9d91582573c' \
  "$patch" | sha256sum -c -

git apply --check "$patch"
git apply "$patch"

echo "Applied the exact frozen-host EEQ block solver on top of persistent EEQ work arrays."
