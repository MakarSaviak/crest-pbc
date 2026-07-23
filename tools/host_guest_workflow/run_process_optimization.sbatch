#!/usr/bin/env bash
#SBATCH --job-name=crest-host-guest-opt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=crest_scheduler_%j.out

set -euo pipefail

: "${CREST_EXE:?set CREST_EXE to the CREST 3.1 executable}"
: "${INPUT_ENSEMBLE:?set INPUT_ENSEMBLE to the input XYZ ensemble}"
: "${CONSTRAINTS:?set CONSTRAINTS to the CREST constraints file}"
: "${OUTPUT_DIR:?set OUTPUT_DIR to the persistent output directory}"
: "${SCRATCH_DIR:?set SCRATCH_DIR to a node-local or project scratch directory}"

workers=${SLURM_CPUS_PER_TASK:-1}
python3 tools/host_guest_workflow/crest_process_scheduler.py \
  --input "$INPUT_ENSEMBLE" \
  --constraints "$CONSTRAINTS" \
  --exe "$CREST_EXE" \
  --output "$OUTPUT_DIR" \
  --scratch-root "$SCRATCH_DIR" \
  --run-id "${SLURM_JOB_ID:-manual}" \
  --workers "$workers" \
  --schedule dynamic \
  --order fifo \
  --frozen-atoms 808 \
  --opt-level vloose \
  --resume
