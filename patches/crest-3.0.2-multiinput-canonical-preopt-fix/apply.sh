#!/usr/bin/env bash
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
patch_dir="$root/patches/crest-3.0.2-multiinput-canonical-preopt-fix"
patch="$patch_dir/multiinput-canonical-preopt-fix.patch"

cd "$root"
printf '%s  %s\n' \
  '7433892e48df721f77b8978e915bc5336409253e5197cffafcf382c60d7caeac' \
  "$patch" | sha256sum -c -

git apply --check "$patch"
git apply "$patch"

echo "Applied CREST 3.0.2 multi-input canonical/preoptimization fix."
