# HISTORICAL — DOES NOT DESCRIBE THE CURRENT SOURCE

This August 2026 implementation note predates capsule V5 and generic
runtime-system validation. It is retained for historical context only; use
`CURRENT_PROCESS_MTD_STATUS.md` for the current supported scope and validation
requirements.

# Dynamic L3-aware process-MTD scheduler — implementation handoff

## Source basis

This tree is derived only from the supplied oracle-fix source archive
`crest-poststage-v2-oraclefix-20260806.tar.gz` (SHA-256
`10492db09e70a673de8cd7dc89538775f91abe32e761f535112814109c018043`).
The optimizer first-touch archive was used as reference only and was not merged.

## Implemented behavior

- Process-MTD scheduling is runtime **N workers x K OpenMP threads/worker**.
- The former process-MTD fixed 48-worker / 4-thread / 192-CPU limits are gone.
- The former GFN-FF `max_inner_threads=4` cap was removed from both MTD entry
  paths. K is derived from CREST's existing `-TMD`/thread budget logic.
- Every supported **multi-trajectory** MTD call routes through the
  process-isolated `crest_search_multimd2` scheduler. The legacy nested-OpenMP
  same-input scheduler is unreachable for `nsim > 1`. A single MTD may retain
  the direct native path.
- The deprecated `CREST_EXPERIMENTAL_MTD_PROCESS_ISOLATION=1` setting is no
  longer required. If the variable is present, only the historical exact value
  `1` is accepted; any other value fails closed.
- Layout validation fails before spawning if N, K, N*K, CPU availability, SMT
  assumptions, or topology discovery are invalid.
- Worker child environment uses runtime K for `OMP_NUM_THREADS` and
  `OMP_THREAD_LIMIT`, explicit singleton `OMP_PLACES`, `OMP_PROC_BIND=close`,
  and BLAS/OpenBLAS/MKL/BLIS threads fixed at 1.

## L3 / NUMA placement

The process-MTD runtime discovers the allowed CPU mask and Linux sysfs topology,
including unified level-3 cache identity from `cache/index*/shared_cpu_list`.
It then:

1. seeds workers across L3 domains first, balancing NUMA nodes and sockets;
2. adds worker lanes round-by-round;
3. prefers the same L3, then the same NUMA node, then the same socket, then any
   remaining allowed core;
4. binds worker memory with `MPOL_BIND` to the NUMA node(s) actually used by
   that worker before capsule deserialization / first touch.

For the measured Fir topology this gives the intended behavior: 48xT1 uses all
24 L3s with two processes/L3; 48xT2 uses all 24 L3s with four active cores/L3;
48xT4 uses all 24 L3s with all eight cores/L3. K up to eight remains inside a
single eight-core L3 when capacity permits; larger K expands within the same
NUMA node before crossing it.

## Oracle / pre-MTD setup

The mandatory production integrity oracle remains fail-closed and is still a
**deterministic T1 inherited-post-trial vs fresh-capsule bitwise oracle**. The
optional fresh-T1 vs fresh-TK numerical reduction-order oracle remains opt-in
only with `CREST_VALIDATE_PROCESS_MTD_THREADS=1`.

The old oracle assumed exactly eight mapped input IDs. It now derives the probe
set from the unique prepared geometries themselves (exact atom order and
bitwise coordinates). Therefore:

- 8 inputs x 6 biases -> 8 unique oracle probes;
- 1 input x 6 biases -> 1 unique oracle probe;
- repeated biases do not cause redundant history probes.

The oracle runs on deep calculator copies, resets the global `engrad_total`
call counter to its entry value, and restores prior OpenMP dynamic/max-thread/
active-level state and the prior OpenBLAS thread count before returning.

## Optimizer scope

Optimizer parallelism is intentionally not generalized. The original optimizer
CPU grouping and gating remain separate from process-MTD L3 discovery. The
post-MTD optimizer-affinity arm is preserved **only** for the historical
validated 48x4, 192-thread full-node case. Other valid N x K MTD batches leave
that optimizer gate disabled.  Once armed, the implicit process-batch gate
enables affinity only for an auto-selected 192x1 optimizer stage; if a later
stage shrinks to another topology (for example 47x4), it falls back to normal
optimizer scheduling without error.  Explicit
`CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY=1` remains strict and rejects any
non-192x1 topology.

## Local validation completed

Using GNU 14 locally and only the bundled GFN-FF + toml-f dependencies:

- full CMake configure and full CREST compile: PASS;
- `mtd_process_runtime.c` `-Wall -Wextra -Werror -fsyntax-only`: PASS;
- capsule-v4 positive/negative round-trip driver with **K=8**: PASS
  (`CAPSULE_NCI_WALLS_EXACT_GATE_PASS`, `CAPSULE_NEGATIVE_GATES_PASS`,
  `CAPSULE_ROUNDTRIP_PASS`);
- process scheduler resolver: absent -> requested=T/status=0; legacy `1` ->
  requested=T/status=0; legacy `0` -> fail-closed configuration status;
- synthetic measured-Fir topology planner oracle: PASS for 48xT1, 48xT2,
  48xT3, 48xT4, 24xT8, 6xT1, and 8xT16; dynamic child env K=8/BLAS=1: PASS;
  over-capacity request: PASS (rejected).

The inherited optimizer-affinity aggregate test runner stops at its existing
caller-source oracle failure in `search_conformers.f90`; running the same runner
against the untouched supplied oracle-fix base fails at the identical source
location. This is therefore not introduced by this scheduler patch.

## Still required on Fir before production use

1. Build this source with the pinned production OpenBLAS and the normal Fir
   compiler/runtime stack.
2. Run a short 1-input x 6-bias AM03@MOF process-MTD job and require the
   mandatory oracle to report 1/1 unique geometries.
3. Run a short 8-input x 6-bias job and require 8/8 oracle probes, 48 clean
   workers, exact trajectories, and the expected L3/NUMA masks.
4. Confirm 48xT1/T2/T4 placement on a Fir node from worker logs; the T1/T2
   patterns should reproduce the successful explicit L3-spread benchmark.
5. Run the known production 8x6 workflow and compare scientific/poststage output
   against the known-good release-v6 reference before replacing the executable.
