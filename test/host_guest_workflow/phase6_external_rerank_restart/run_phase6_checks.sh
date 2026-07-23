#!/usr/bin/env bash
set -euo pipefail
exe=${1:?usage: run_phase6_checks.sh /path/to/crest}
root=$(cd "$(dirname "$0")" && pwd)
work="$root/run"
rm -rf "$work"
mkdir -p "$work/base"
cp "$root/input.toml" "$root/struc.xyz" "$work/base/"

# Stage 1: execute one of the requested two iterations and pause deliberately.
(
  cd "$work/base"
  "$exe" input.toml -nomtd > first.out 2> first.err
)
grep -F 'Meta-Dynamics Iteration 1' "$work/base/first.out"
grep -F 'External-rerank checkpoint reached after iteration 1.' "$work/base/first.out"
grep -F 'CREST terminated normally.' "$work/base/first.out"
grep -F 'stage awaiting_external_rerank' "$work/base/crest.restart"
grep -F 'mtd_iter 1' "$work/base/crest.restart"
grep -F 'target_mtd_iter 2' "$work/base/crest.restart"
grep -F 'settings_file input.toml' "$work/base/crest.restart"
[[ -s "$work/base/.cre_0.xyz" ]]

# Stage 2: the positional XYZ supplies only the replacement geometry. The TOML
# is restored from the checkpoint, including GFN-FF fragments.
cp -a "$work/base" "$work/valid"
cp "$work/valid/.cre_0.xyz" "$work/valid/crest-best-external.xyz"
(
  cd "$work/valid"
  "$exe" crest-best-external.xyz --restart > second.out 2> second.err
)
grep -F 'External-rerank seed validation:' "$work/valid/second.out"
grep -F 'GFN-FF fragment groups used: 2' "$work/valid/second.out"
grep -F 'inter-fragment contacts were excluded from the topology comparison' "$work/valid/second.out"
grep -F 'coordinate equality was not required' "$work/valid/second.out"
grep -F '# fragment in coord' "$work/valid/second.out"
! grep -q 'Meta-Dynamics Iteration 1' "$work/valid/second.out"
grep -F 'Meta-Dynamics Iteration 2' "$work/valid/second.out"
grep -F 'CREST terminated normally.' "$work/valid/second.out"
grep -F 'stage done' "$work/valid/crest.restart"

# Break one O-H bond inside fragment 1 while retaining atom count/order and the
# extxyz schema. Validation must fail before iteration 2 begins.
cp -a "$work/base" "$work/broken_topology"
python3 - "$work/broken_topology/.cre_0.xyz" \
  "$work/broken_topology/crest-best-external.xyz" <<'PY'
from pathlib import Path
import sys
src, dst = map(Path, sys.argv[1:])
lines = src.read_text().splitlines()
n = int(lines[0])
out = lines[:2]
for atom_index, line in enumerate(lines[2:2+n], start=1):
    fields = line.split()
    values = list(map(float, fields[1:]))
    if atom_index == 2:
        values[0] += 4.0
    out.append(f"{fields[0]:<2s} " + " ".join(f"{v:20.12f}" for v in values))
out.extend(lines[2+n:])
dst.write_text("\n".join(out) + "\n")
PY
set +e
(
  cd "$work/broken_topology"
  "$exe" crest-best-external.xyz --restart > second.out 2> second.err
)
rc=$?
set -e
[[ $rc -ne 0 ]]
! grep -q 'Meta-Dynamics Iteration 2' "$work/broken_topology/second.out"
grep -F 'external-rerank seed topology does not match the iteration-1 ensemble' \
  "$work/broken_topology/second.err"

echo 'phase6 external-rerank restart checks passed'
