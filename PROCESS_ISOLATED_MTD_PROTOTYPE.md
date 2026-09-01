# HISTORICAL — DOES NOT DESCRIBE THE CURRENT SOURCE

This archived V9 prototype description predates the current generic
implementation. Its fixed 938-atom/48-job contract, capsule-V4 description,
and opt-in scheduler wording are historical only. Do not use it for a release
decision; see `CURRENT_PROCESS_MTD_STATUS.md`.

# Release-v6 combined fragment-box and optimizer-affinity prototype V9

Status: local source prototype only. It has not been deployed, submitted,
scientifically accepted, or promoted.

## Immutable starting point

This V9 tree is a documentation-only successor to the source-frozen V8
combined prototype:

```text
/scratch/saviak/crest-eeq-fir-secondary-20260729/prototypes/release_v6_process_isolated_mtd_scheduler_always_real64_eeq_fitde_isolated_v6_statusprop_nciwalls_fragmentbox_optspread_v8_20260730
```

The complete V8 input tree SHA-256 was
`b007fb57354cb5d47b4db9a53e95bad418f669e0d7d19bad64dc4fbe35c1c226`.
V9 changes documentation only; every production and prototype-test source
file is byte-identical to V8. V8 itself combined the frozen V7 affinity source
and boxed-fragment manifest
`4e0a6b7521abc65847107c3a515d45c1cf861dd37ee993679a32f29bf58aae91`.
The exact composition, decomposition invariants, and remaining gates are in
`COMBINED_CANDIDATE_V9.md`. The rest of this document describes the inherited
V7 process scheduler and affinity behavior; V9 additionally contains only the
fragment-string representation and its serialization/access sites.

## Exact integration seam

The primary production-flow dispatch is inside `crest_search_multimd2()` in
`src/algos/parallel.f90`, after the parent has already prepared the mapped
`mols(1:48)` and resolved `mddats(1:48)`. The process path replaces only that
routine's OpenMP dynamics-execution block. Additional guards reject the opt-in
from the legacy same-input scheduler and hard-stop on invalid configuration or
process failure. A common post-branch status gate in the conformer-search
caller also returns before any trajectory consumer when an MTD scheduler sets
failure status. Together these prevent native fallback or stale-trajectory
optimization.

Optimizer affinity additionally edits the surrounding control flow of
`crest_oloop()` in `src/algos/parallel.f90`: it resolves the gate, prepares and
binds the outer team, barriers before task creation, restores the parent mask,
and returns typed failure status after cleanup. The 60-line scientific OpenMP
task body itself remains byte-for-byte identical to the release-v6 starting
source. V7 also adds immediate typed-status return guards to every direct
`crest_oloop()` caller and every caller that can receive a nested
`crest_refine()` failure. Those guards run before caller output, sorting,
renaming, or other postprocessing; the normal-success path is unchanged.

The parent still owns the following scientific behavior unchanged:

- all raw-input/TOML parsing and eight-input x six-bias mapping;
- trial MTD and every nonmatching MD/MTD mode;
- the release-v6 numeric-order, exact-byte trajectory collector;
- the 60-line optimization task body and all optimization kernels;
- success-path screening, CREGEN, and final output;
- external-rerank checkpoint/return/restart behavior.

The default is byte-for-byte release-v6 control flow when the environment
variable is absent. The process path is entered only for the exact value
`CREST_EXPERIMENTAL_MTD_PROCESS_ISOLATION=1`. Any other present value,
including an overlong value, fails closed; it never silently selects the
native nested scheduler. Legacy same-input/later-iteration MTD entry points
also reject the opt-in, so a mispackaged one-input target cannot fall back to
native scheduling.

The v1 scope is deliberately narrow:

- exactly 48 jobs and 192 physical CPUs;
- exactly four built-in GFN-FF threads per trajectory;
- exactly 938 atoms, with atoms 1--808 frozen and atoms 809--938 mobile;
- exactly one active GFN-FF calculator level;
- exactly two GFN-FF fragments, `1-808` and `809-938`;
- exactly one occurrence of every input/bias pair in the 8 x 6 map;
- already-resolved stochastic RMSD MTD with a positive mass-weighted COM bias
  and an exact guest inclusion mask for atoms 809--938;
- `shake=false`, no prepared SHAKE state, and no MD restart/history;
- an explicit existing canonical topology restart;
- no `samerand`, mixed/external calculators, ONIOM, calculator constraints,
  scans, numerical gradients, or nonempty PDB payloads.

## Capsule semantics

`src/mtd_process_capsule.F90` writes the already-prepared molecule, resolved
`mddata`/`mtdpot` state, freeze mask, fragments, and restricted GFN-FF
calculator settings to one capsule per worker. Workers do not parse or reread
the original coordinate ensemble or TOML. At opt-in batch entry, the parent
calls the same idempotent `mddata%defaults()` routine that native `dynamics()`
calls before use. This defines upstream `maxblock` and resolves any remaining
fallback values before serialization; the child repeats the same operation.

The exact target has `shake=false`. Active/initialized SHAKE, constraint/WBO
arrays, per-step SHAKE scratch, and a live SHAKE freeze pointer are rejected.
The upstream SHAKE pointer lacks a default null initializer, so the opt-in
entry explicitly nullifies that semantically dead pointer before testing or
copying it. Undefined SHAKE `dro`/`dr` scratch bytes are therefore never read
or serialized.

Capsule v4 is a native-endian binary format for parent and child instances of
the same executable. It has a fixed magic, version, endian sentinel, bounded
array dimensions, explicit integer widths, logical-byte encoding, and a fixed
trailer with no permitted trailing data. `wp` values are transferred as their
native real64 bytes; there is no decimal conversion or precision loss. It is
not an archival or cross-compiler interchange format.

Before any child starts, every capsule is decoded, its full supported state is
compared against the parent state, re-encoded, and required to be
byte-identical. The reader rejects truncated, wrong-endian, malformed-shape,
and trailing-byte capsules before dynamics.
`prototype_tests/capsule_roundtrip_driver.F90` supplies a separate 938-atom,
SHAKE-disabled, COM-biased/frozen-host sentinel and negative cases for these
contracts.

Capsule v4 intentionally does not serialize `gfnff_data`: it rejects any
calculator settings object whose `ff_dat` runtime state is allocated. Normal
release-v6 trial MTD does leave that state allocated, so the scheduler handles
the boundary explicitly before writing a capsule. It deep-copies the parent
calculator, verifies the inherited restart provenance, initialized topology
and neighbour list, exact fragment and freeze masks, and absence of the
nonserialized `sTorsl` term, then deallocates `ff_dat` on that independent copy
only. It also restricts the one-level selector to weighted (`id=0`) or direct
(`id=1`) mode. In weighted mode, any inherited one-element `eweight` cache must
be bitwise equal to the configured level weight. The clean copy drops
calculator energy/gradient/weight and backup scratch so its first call
reconstructs exactly the state that the capsule decoder creates. The live
parent calculator is never changed or deallocated.

That clean copy is accepted only after a fail-closed, pre-spawn calculator
oracle. For each of the eight input structures, an independent clone of the
untouched post-trial calculator is compared with an independent clone freshly
reconstructed from the byte-preserved topology. Both make a four-thread first
energy/gradient call and a second call after an exactly representable
`2^-7`-Bohr displacement of guest atom 809. The oracle requires bitwise-equal
energy, all 20 reported GFN-FF result values, all charges, the full gradient,
and the complete mutable HB/XB neighbour-list flags, counts, reference
geometry, and every active list entry. GFN-FF reads `hblist1`, `hblist2`, and
`hblist3` only through `nhb1`, `nhb2`, and `nxb`; their larger allocated second
dimensions are spare capacity and can contain stale, unread tails after a
rebuild. The oracle therefore requires valid capacity on both sides but does
not mistake capacity size or inactive tail bytes for scientific state. It also
requires the same topology-restart provenance.
The 32 probe calls are removed from CREST's global call counter, and the
parent's prior OpenMP dynamic/thread/active-level and OpenBLAS thread settings
are restored. Any mismatch stops the opt-in path before a worker is spawned;
there is no tolerance, fallback, or silent state loss.

## Process, topology, and output isolation

`src/mtd_process_runtime.c` invokes the same executable using `posix_spawn()`;
the internal worker mode is selected before normal CREST argument parsing.
There is no `fork()` followed by execution of Fortran/OpenMP code.

All artifacts live under the persistent private directory
`crest_process_mtd_audit/process_batch_<parent-pid>`, retained on success or
failure. The audit root, batch directory, worker directory, and calculator
directory are verified as owner-only mode `0700`. Every worker receives a
private directory containing:

- its capsule and result record;
- a byte copy of the parent's explicitly configured canonical `gfnff_topo`;
- its private MD restart and trajectory;
- a private calculator directory, stdout, and stderr;
- any other relative files emitted by GFN-FF in that private working directory.

The parent first makes and byte-verifies one preserved topology master, then
makes one byte-verified private copy per worker. The worker sets its calculator
restart path to that private copy, and the parent verifies after completion
that every worker copy is still byte-identical to the master. No GFN-FF, EEQ,
frozen-host validation/cache, HB/XB, fragment, COM-bias, or dynamics scientific
kernel is edited. Each process constructs private topology/workspace state
from that byte-identical restart only after the pre-spawn oracle has proved the
fresh state bitwise equivalent for the exact eight-input target; no topology
or writable calculator object is shared between workers.

## Affinity and threading

The launcher requires exactly 192 online and allowed physical CPUs, a Slurm
full-node allocation marker, and an affinity mask containing exactly those
192 CPUs. It reads package, die, core, and NUMA identity for every allowed CPU
from `/sys`, fails closed if any topology field is unavailable, and rejects
SMT siblings masquerading as independent cores. CPUs are grouped four at a
time within one NUMA node, and the 48 groups are ordered round-robin across
NUMA nodes: spread between workers, close within each worker.

CPU topology is scanned and cached once by the parent. Before `exec`, each
child environment receives an explicit four-CPU
`OMP_PLACES`, `OMP_PROC_BIND=close`, `OMP_NUM_THREADS=4`, one active OpenMP
level, and `OPENBLAS_NUM_THREADS=1` (plus equivalent BLAS variables). Inherited
`GOMP_CPU_AFFINITY`, `GOMP_SPINCOUNT`, `OMP_SCHEDULE`, `OMP_WAIT_POLICY`, and
conflicting OpenMP/BLAS values are removed. Before decoding or first-touching
the scientific capsule, the child narrows its Linux affinity mask to the same
four CPUs and binds memory allocation to their NUMA node. It then runs a
four-place/four-thread self-check that verifies every actual CPU and singleton
OpenMP place.

Production RNG remains stochastic. The parent refuses any resolved MTD with
`samerand=true`; the prototype does not add a deterministic seed.

## Status, validation, and cleanup

The parent records both the real `waitpid()` exit/signal status and the
`dynamics()` termination status returned in a versioned, exact-size private
result file. The result also carries that child's energy/gradient-call count;
the parent accepts only exactly `length_steps` calls per child and adds the
validated aggregate to its own counter.

A successful child is not released to collection until both logs contain
exactly one normal-termination marker in total and no CREST, native-signal,
Fortran, allocator, BLAS, loader, solver, IEEE, or NaN fatal signature. Its
trajectory must have the exact expected frame count, complete 938-atom XYZ
records in the original atom-symbol order, finite `Epot` and coordinates, and
nonzero size. Only after all 48 workers pass does the unchanged release-v6
numeric-order byte-stream collector run. The parent emits native-compatible
`*MTD ... completed successfully` markers for the persistent auditor.

On spawn or child failure, the parent sends SIGTERM to every unreaped child,
then uses bounded reaping with SIGKILL fallback rather than waiting forever.
Each child sets Linux `PR_SET_PDEATHSIG=SIGKILL` and verifies its expected
parent. SIGKILL is deliberate because CREST replaces SIGTERM with a Fortran
handler during dynamics; abrupt parent death must not depend on that handler
or active OpenMP/allocator state.
Failure artifacts and private directories are preserved for audit. Invalid
opt-in state or a process-batch failure hard-stops at the scheduler boundary;
the common conformer-search status gate additionally blocks every returned
MTD failure before downstream collection or optimization. Optimizer-affinity
setup, bind, restore, or finish failure returns from `crest_oloop()` with typed
status, and every audited direct or transitive caller immediately returns
before caller output or postprocessing. A stale trajectory or failed
optimization result cannot be consumed.

## Changed files

- `src/mtd_process_runtime.c` (new)
- `src/mtd_process_capsule.F90` (new)
- `src/mtd_process_scheduler.F90` (new)
- `src/algos/parallel.f90`
- `src/algos/search_conformers.f90`
- `src/algos/search_1.f90`
- `src/algos/optimization.f90`
- `src/algos/refine.f90`
- `src/algos/search_mecp.f90`
- `src/algos/protonate.f90`
- `src/calculator/calc_type.f90`
- `src/calculator/api_helpers.F90`
- `src/parsing/parse_calcdata.f90`
- `src/external_rerank_restart.f90`
- `src/crest_main.f90`
- `src/CMakeLists.txt`
- `src/meson.build`
- `prototype_tests/capsule_roundtrip_driver.F90` (new)
- `prototype_tests/process_option_driver.F90` (new)
- `prototype_tests/optimizer_affinity_c_driver.c` (new)
- `prototype_tests/optimizer_affinity_fortran_driver.F90` (new)
- `prototype_tests/check_optimizer_affinity_callers.py` (new in V7)
- `prototype_tests/run_optimizer_affinity_oracles.sh` (new)
- `prototype_tests/fragment_string_copy_oracle.F90` (new in V8)
- `prototype_tests/fragment_string_real_copy_oracle.F90` (new in V8)
- `prototype_tests/compute_source_tree_sha256_v1.sh` (pinned provenance helper)
- prototype documentation/evidence files

No GFN-FF, EEQ, HB/XB, fragment-topology, freeze-validation, RMSD, COM-bias,
dynamics, collector, or external-rerank scientific kernel was edited. V8
changes the owning representation and exact serialization/access of fragment
strings because GCC 12 miscompiled the predecessor's nested deferred-length
character array. `crest_oloop()` control flow and caller failure propagation
were edited; its exact 60-line scientific optimization task body remains
byte-identical.

## Checks inherited from V8 and rerun for V9

- exact-source caller oracle inventory for all 14 direct `crest_oloop()` sites
  and all five transitive `crest_refine()` propagation sites;
- regular and AddressSanitizer/UndefinedBehaviorSanitizer C lifecycle oracles;
- regular and sanitizer Fortran control-flow oracles, including zero tasks on
  bind/team failure and no task or output continuation after returned failure;
- production-object proof that test-only injection symbols are absent;
- C11 `-Wall -Wextra -Werror` runtime compilation and gfortran syntax/object
  compilation of every V7-modified Fortran caller against the preserved V6
  module set;
- byte-identity checks for the 60-line optimization task body and inherited
  GFN-FF tree, plus a V6-preservation hash check;
- a complete local affinity-oracle evidence manifest; the oracle evidence
  hashes every scanned or compiled local source input and preserves separate
  stdout and stderr for each check.

These checks are compile/control-flow oriented only. No full build, SSH,
Slurm submission, compute-node run, or scientific acceptance was performed.

## Unresolved blockers before any deployment

1. A clean isolated build from the V9 source manifest with a fully pinned
   compiler, OpenBLAS-256, and linker is still required.
2. The capsule and opt-in sentinels must be linked to that exact clean build
   and rerun; local link checks do not substitute for the clean-build evidence.
3. The 192-core worker affinity/self-check must be exercised on a compute node;
   login-node topology is intentionally unsupported.
4. The real post-trial state must pass the integrated eight-input/two-probe
   bitwise calculator oracle before any child spawn. A short 8 x 6 smoke must
   then prove 48 clean exits, private topology identity, exact call/frame
   counts, no fatal signatures, complete finite trajectories,
   native-compatible completion markers, and a clean unchanged release-v6
   collector.
5. Fixed-topology T1/T2/T4 energy/component/charge/full-gradient correctness
   must pass to scientifically justified reduction-order roundoff. No source
   tolerance or calculator approximation may be added.
6. Order-balanced stochastic short-MTD timing must remeasure integrated
   one-, two-, and four-thread scaling. The current final-like reference is
   approximately 1.63x for T4/T1; the older approximately 2.9x observation was
   not reproduced and is not acceptance evidence for this candidate.
7. Fixed 54,300-frame screen and one-iteration production-like gates remain
   mandatory before promotion.
8. The prototype deliberately reloads each private topology from the exact
   configured canonical restart file; it does not serialize arbitrary mutated
   in-memory `gfnff_data`. Capsule v4 rejects allocated `ff_dat`; the scheduler
   emits clean settings only after the exact post-trial state passes its
   bitwise oracle. Therefore v1 rejects configurations without an explicit
   topology restart and must be tested only against the pinned canonical
   topology contract.
9. The prototype is Linux/Slurm specific: it requires a full 192-physical-core
   node, `/sys` CPU/NUMA topology, `sched_setaffinity`, Linux memory policy, and
   `PR_SET_PDEATHSIG`. Failure of any binding proof is a hard rejection.
