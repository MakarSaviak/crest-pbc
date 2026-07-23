#!/usr/bin/env bash
set -euo pipefail
exe=${1:?usage: run_phase1_checks.sh /path/to/crest}
root=$(cd "$(dirname "$0")" && pwd)
rm -rf "$root/run_cli" "$root/run_toml" "$root/run_default" "$root/run_invalid"
mkdir -p "$root/run_cli" "$root/run_toml" "$root/run_default" "$root/run_invalid"

(
  cd "$root/run_cli"
  cp "$root/struc.xyz" "$root/input_cli.toml" .
  "$exe" input_cli.toml -TMD 2 > crest.out 2>&1
)

grep -F "MTD/MD runs: 2 parallel jobs × 1 core/job  (2 threads" "$root/run_cli/crest.out"
grep -E "optimizations: .*\(4 threads" "$root/run_cli/crest.out"

(
  cd "$root/run_toml"
  cp "$root/struc.xyz" "$root/input_toml.toml" .
  "$exe" input_toml.toml > crest.out 2>&1
)

grep -F "MTD/MD runs: 2 parallel jobs × 1 core/job  (2 threads" "$root/run_toml/crest.out"
grep -E "optimizations: .*\(4 threads" "$root/run_toml/crest.out"


(
  cd "$root/run_default"
  cp "$root/struc.xyz" "$root/input_cli.toml" .
  "$exe" input_cli.toml > crest.out 2>&1
)

grep -F "MTD/MD runs: 4 parallel jobs × 1 core/job  (4 threads" "$root/run_default/crest.out"
grep -E "optimizations: .*\(4 threads" "$root/run_default/crest.out"

set +e
(
  cd "$root/run_invalid"
  cp "$root/struc.xyz" "$root/input_cli.toml" .
  "$exe" input_cli.toml -dry -TMD 0 > zero.out 2>&1
); rc0=$?
(
  cd "$root/run_invalid"
  "$exe" input_cli.toml -dry -TMD 1.5 > fractional.out 2>&1
); rcf=$?
(
  cd "$root/run_invalid"
  "$exe" input_cli.toml -dry -TMD > missing.out 2>&1
); rcm=$?
set -e

[[ $rc0 -ne 0 && $rcf -ne 0 && $rcm -ne 0 ]]
grep -F -- "-TMD requires a positive integer" "$root/run_invalid/zero.out"
grep -F -- "-TMD requires a positive integer" "$root/run_invalid/fractional.out"
grep -E "missing|argument" "$root/run_invalid/missing.out"

echo "phase1 -TMD checks passed"
