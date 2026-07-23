# Process-isolated CREST 3.1 ensemble optimization

`crest_process_scheduler.py` runs each input frame in a fresh, single-threaded
CREST process and private scratch directory.  It was adapted to the CREST 3.1
CLI (`crest --mdopt input.xyz ...`) and is intended for large frozen-host
host–guest ensembles where optimizer state, calculator state, and temporary
files must never be shared between frames.

## Example

```bash
python3 tools/host_guest_workflow/crest_process_scheduler.py \
  --input ensemble.xyz \
  --constraints constraints.inp \
  --exe /path/to/crest \
  --output results \
  --scratch-root /scratch/$USER/crest \
  --run-id am03-mof5 \
  --workers 32 \
  --schedule dynamic \
  --order fifo \
  --frozen-atoms 808 \
  --opt-level vloose \
  --overwrite
```

Every worker process receives one OpenMP/BLAS thread.  Validated final
structures are reconstructed in the original input-frame order regardless of
completion order.

## Scheduling

- `dynamic`: workers draw frames from a common queue.
- `static-contiguous`: each worker receives a contiguous frame block.
- `static-roundrobin`: frames are distributed cyclically.
- `--order predicted-longest-first --predicted-costs FILE`: launch expensive
  frames first while preserving original ordering in the final ensemble.

## Resume and validation

`--resume` reuses a frame only when its input digest, output digest, convergence
record, and manifest still agree.  A clean one-worker result directory can be
supplied through `--reference` to compare call count, energy, guest RMSD, and
optionally byte identity.

Persistent output contains per-frame JSON/XYZ/log files, `frames.csv`,
`frames.json`, `summary.json`, and the reconstructed `crest_ensemble.xyz`.
Transient CREST workspaces are removed after validation unless
`--keep-workspaces` is set.

## CREST 3.1 frozen-host benchmark

For four representative 938-atom AM03@MOF frames with atoms 1–808 frozen, the
reduced optimizer completed 254 total energy-gradient calls in 62.95 s with one
worker and 25.82 s with four workers.  Results matched clean CREST 3.1 within
`1.0e-10 Eh` and `7.9e-9 Å`, with zero host displacement.
