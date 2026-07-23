#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cat > "$work/fake_crest.py" <<'PY'
#!/usr/bin/env python3
from pathlib import Path
src=Path('input.xyz').read_text().splitlines()
# Preserve coordinates while normalizing the energy comment into a CREST-like one.
src[1]='-1.000000000000'
Path('crest_ensemble.xyz').write_text('\n'.join(src)+'\n')
print('1 of 1 structures successfully optimized')
print('Total number of energy+grad calls: 3')
PY
chmod +x "$work/fake_crest.py"
cat > "$work/in.xyz" <<'XYZ'
2
frame=0
H 0 0 0
H 0 0 1
2
frame=1
H 0 0 0
H 0 0 2
2
frame=2
H 0 0 0
H 0 0 3
2
frame=3
H 0 0 0
H 0 0 4
XYZ
cat > "$work/constraints.inp" <<'EOF2'
$fix
 atoms: 1
$end
EOF2
python3 "$root/crest_process_scheduler.py" \
  --input "$work/in.xyz" --constraints "$work/constraints.inp" \
  --exe "$work/fake_crest.py" --output "$work/out" \
  --scratch-root "$work/scratch" --workers 2 --frozen-atoms 1 --overwrite
python3 - "$work/out" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])
s=json.loads((p/'summary.json').read_text())
f=json.loads((p/'frames.json').read_text())
assert s['completed_frames']==4 and s['validation_passed']
assert s['total_calls']==12
assert [x['frame'] for x in f]==[0,1,2,3]
assert all(x['host_max_displacement_a']==0 for x in f)
assert (p/'crest_ensemble.xyz').read_text().count('\n2\n')==3
PY
python3 "$root/crest_process_scheduler.py" \
  --input "$work/in.xyz" --constraints "$work/constraints.inp" \
  --exe "$work/fake_crest.py" --output "$work/out" \
  --scratch-root "$work/scratch" --workers 2 --frozen-atoms 1 --resume
python3 - "$work/out/summary.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
assert s['resumed_frames']==4 and s['completed_frames']==4
assert s['timing_status']=='resume_segment_only'
PY
echo 'process scheduler smoke test passed'
