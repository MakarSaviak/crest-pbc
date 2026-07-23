#!/usr/bin/env bash
set -euo pipefail
exe=${1:?usage: run_phase3_checks.sh /path/to/crest}
root=$(cd "$(dirname "$0")" && pwd)
run_common() {
  local input=$1 outdir=$2
  rm -rf "$outdir"
  mkdir -p "$outdir"
  cp "$root/$input" "$root/mtdbias" "$outdir/"
  (
    cd "$outdir"
    "$exe" "$input" --gfnff --imtdgc -nci -readbias -mrest 1 \
      -len 0.005 -mddump 1 -T 2 -TMD 2 > crest.out 2>&1
  )
}

run_common placements.xyz "$root/run_multi"
grep -F "Detected 2 positional XYZ frames" "$root/run_multi/crest.out"
grep -F "NCI first MTD batch: 2 inputs x 2 biases = 4 trajectories" "$root/run_multi/crest.out"
grep -F "CREST terminated normally" "$root/run_multi/crest.out"
[[ -s "$root/run_multi/crest_nci_input_starts.xyz" ]]
[[ -s "$root/run_multi/crest_nci_mtd_jobs.tsv" ]]
python3 - "$root/run_multi/crest_nci_mtd_jobs.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    rows=list(csv.DictReader(f, delimiter='\t'))
assert len(rows)==4, rows
expected=[('1','1'),('1','2'),('2','1'),('2','2')]
assert [(r['input_placement'],r['bias_configuration']) for r in rows]==expected
assert [r['job_index'] for r in rows]==['1','2','3','4']
assert all(r['termination_status']=='0' for r in rows)
assert len({r['trajectory'] for r in rows})==4
PY

run_common single.xyz "$root/run_single"
! grep -q "Detected .* positional XYZ frames" "$root/run_single/crest.out"
! grep -q "NCI first MTD batch" "$root/run_single/crest.out"
[[ ! -e "$root/run_single/crest_nci_mtd_jobs.tsv" ]]
grep -F "CREST terminated normally" "$root/run_single/crest.out"
echo "phase3 multi-input NCI checks passed"
