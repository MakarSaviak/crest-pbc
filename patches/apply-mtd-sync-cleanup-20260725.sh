#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch="$root/patches/mtd-sync-cleanup-20260725.patch"

cd "$root"
printf '%s  %s\n' \
  '16b6f56c0d7f7550ee1ba8bb626c1ba677c7d7b1b927c1272257401aef145550' \
  "$patch" | sha256sum -c -

git apply --check "$patch"
git apply "$patch"

echo "Applied independent-trajectory MTD synchronization cleanup on top of the exact EEQ solver."
