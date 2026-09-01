# Current generic process-isolated MTD status

This document describes the current source tree. It supersedes the historical
V8/V9 process-MTD notes retained alongside it for provenance.

## Supported scope

The process scheduler is system-generic for its supported workflow: a
multi-trajectory GFN-FF NCI/RMSD-MTD run. It obtains the atom count from the
runtime molecule; accepts arbitrary atom order, freeze masks, MTD inclusion
masks, GFN-FF fragment strings/counts, and one automatic all-atom NCI wall;
and uses capsule V5. It is not a general serialization layer for every CREST
MTD feature. In particular, the process route deliberately rejects active
SHAKE/restarts and unsupported calculator or MTD state rather than dropping it.

## Capsule and oracle safety

Capsule V5 has `CREST_MTD_PROCESS_CAPSULE_V5` magic and rejects earlier-format
headers. Fragment strings retain their individual allocated lengths. Empty or
zero-length allocated fragment specifications fail closed; this is an
intentional validation improvement over the former implicit behavior.

The parent calculator oracle selects a movable included atom where possible;
it never deliberately perturbs an atom marked frozen. If every atom is frozen,
the controlled-motion portion is skipped with a diagnostic.

## Validation-only work

The V4-header regression and prototype drivers are build/CI validation only.
They are not invoked by normal production MTD jobs.

## Release gates

Before promotion, run a clean OpenBLAS build reporting `MAX_THREADS >= 256`,
confirm children run with `OPENBLAS_NUM_THREADS=1`, and execute real small
non-AM03 and AM03 multi-trajectory process-MTD regressions. A source archive
for publication must exclude compiled binaries and build/evidence directories.
