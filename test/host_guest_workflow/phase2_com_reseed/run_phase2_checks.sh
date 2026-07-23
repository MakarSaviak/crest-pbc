#!/usr/bin/env bash
set -euo pipefail
exe=${1:?usage: run_phase2_checks.sh /path/to/crest}
root=$(cd "$(dirname "$0")" && pwd)
rm -rf "$root/run_toml" "$root/run_xtb" "$root/run_invalid"
mkdir -p "$root/run_toml" "$root/run_xtb" "$root/run_invalid"
(
  cd "$root/run_toml"
  cp "$root/struc.xyz" "$root/input_toml.toml" .
  "$exe" input_toml.toml > crest.out 2>&1
)
grep -E "COM bias[[:space:]]*:.*enabled" "$root/run_toml/crest.out"
grep -E "COM factor /[[:space:]]*Eh" "$root/run_toml/crest.out"
grep -F "COM mass weighted" "$root/run_toml/crest.out"
(
  cd "$root/run_xtb"
  cp "$root/struc.xyz" "$root/input_xtb.toml" "$root/metadyn.inp" .
  "$exe" input_xtb.toml -cinp metadyn.inp -TMD 1 > crest.out 2>&1
)
grep -E "COM bias[[:space:]]*:.*enabled" "$root/run_xtb/crest.out"
grep -E "COM factor /[[:space:]]*Eh" "$root/run_xtb/crest.out"
grep -F "COM mass weighted" "$root/run_xtb/crest.out"
set +e
(
  cd "$root/run_invalid"
  cp "$root/struc.xyz" "$root/input_invalid.toml" .
  "$exe" input_invalid.toml > crest.out 2>&1
); rc=$?
set -e
[[ $rc -ne 0 ]]
grep -F "COM MTD width must be positive" "$root/run_invalid/crest.out"
# The two production search drivers must refresh mol from env%ref immediately
# before constructing each new MTD batch.
grep -A5 -F "Always seed a new MTD iteration" "$root/../../../src/algos/search_conformers.f90" | grep -F "call env%ref%to(mol)"
grep -A5 -F "Always seed a new MTD iteration" "$root/../../../src/algos/search_entropy.f90" | grep -F "call env%ref%to(mol)"
echo "phase2 COM/reseed checks passed"
