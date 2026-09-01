# Local prototype tests

`capsule_roundtrip_driver.F90` keeps the 938-atom AM03-shaped input as a
regression fixture, then exercises two generic cases: a 17-atom no-freeze case
with COM bias disabled, and a 47-atom case with a non-contiguous freeze mask,
three independently sized fragment strings (including one longer than 16
characters), and an inclusion mask unrelated to the frozen atoms. Every case
uses SHAKE-disabled stochastic RMSD-MTD, an explicit topology restart, and one
marked automatic NCI wall sized from its molecule. It writes capsule v5,
decodes it, compares supported state, rewrites it, requires byte identity, and
checks molecule, MTD/COM, fragment, freeze-mask, and wall fields.

The same driver requires rejection of:

- a one-byte truncation;
- a damaged endian sentinel;
- a V4 magic/version header presented to the V5 reader;
- one byte appended after the canonical trailer;
- a source-versus-decoded metadata mismatch;
- preallocated dynamics block scratch;
- active SHAKE in this `shake=false` target;
- a zero trajectory-dump interval and nonfresh MTD counter;
- an empty or zero-length fragment array, or an empty/unallocated fragment
  string;
- unsupported reference-path state;
- an allocated GFN-FF `ff_dat` runtime object, which capsule v5 must never
  serialize incompletely or silently discard;
- an allocated calculator weight cache, which the scheduler must validate and
  remove only from its independent clean copy before capsule serialization.

Link it to the prototype's freshly built static CREST library plus the exact
prototype capsule module; do not link it against a different source snapshot.
It prints `CAPSULE_NEGATIVE_GATES_PASS` and `CAPSULE_ROUNDTRIP_PASS` only after
all checks pass.

`process_option_driver.F90` exercises the process-scheduler environment resolver.
The process scheduler is now the default for every supported multi-MTD batch, so
with the deprecated variable absent it must report `requested=T,status=0`; exact
legacy value `1` is still accepted for compatibility. Any other present value
(including `0` and a value longer than 32 bytes) must return nonzero
configuration status. This proves an invalid legacy setting cannot silently
fall back to native scheduling.

The end-to-end smoke is intentionally not a login-node test. On a compute node,
run the reviewed AM03@MOF multi-MTD harness without requiring any process-MTD
opt-in. Before any child is spawned, the scheduler must report that every
unique prepared input geometry passes the mandatory two-step, deterministic-T1
post-trial-inherited versus fresh-restart bitwise oracle. For the established
8-input x 6-bias batch this is 8/8; for a 1-input x 6-bias batch it is 1/1.
This compares energy, all reported
components, charges, the full gradient, and mutable HB/XB neighbour-list state
on each first input frame and an exact `2^-7`-Bohr displacement selected from
the prepared MTD inclusion mask or freeze state. HB/XB
comparison includes all live list entries while excluding only unused
over-allocation tails that the GFN-FF kernels never read. The
scheduler uses the prepared GFN-FF state, must not set `samerand`, and
currently fails closed for active SHAKE because that state is not serialized.
It must then report all N exit-zero workers with exact call and trajectory-frame
counts, clean logs, immutable private topologies, and finite ordered
runtime-sized trajectories before the unchanged collector runs. Worker count N
and threads per worker K are runtime values rather than fixed 48 x 4 constants.

No job script is included here: this source prototype is local-only and must
first pass a clean isolated build, source-manifest pinning, and independent
scientific review before deployment or submission is considered.


## Dynamic L3 planner oracle

`mtd_l3_planner_driver.c` directly exercises the production planner against a
synthetic copy of the measured Fir topology: 192 physical cores, 24 eight-core
L3 domains, eight NUMA nodes, and two sockets. It requires L3-first spreading
for 48xT1/T2/T3/T4, one process per L3 for 24xT8, broad spreading for 6xT1,
and same-NUMA expansion for T16. It also verifies that K=8 is propagated into
the worker OpenMP environment while BLAS remains at one thread, and that an
over-capacity N x K request fails. Compile from the source root with:

```sh
gcc -std=gnu11 -O0 -Wall -Wextra -Wno-unused-function \
  prototype_tests/mtd_l3_planner_driver.c -o /tmp/mtd_l3_planner_driver
/tmp/mtd_l3_planner_driver
```

## Optimizer-affinity local oracles

`optimizer_affinity_c_driver.c` exercises the C lifecycle directly with four
real pthreads: normal bind/restore, active-plan re-entry rejection, one
injected bind failure, a reported-team mismatch, one injected restore failure,
and a second successful region after all failure cases.  The injection entry
point exists only when `mtd_process_runtime.c` is compiled with
`CREST_OPTIMIZER_AFFINITY_TESTING`; the runner verifies that the symbol is
absent from a production object.

The same C driver exercises the production gate.  With no opt-in and no
validated process batch, optimizer behavior stays unchanged.  The exact
standalone-screen opt-in is
`CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY=1`; any other present value fails
closed, and Fortran permits that request only when `env%crestver` is
`crest_screen`.  The real process-child environment constructor removes this
parent-only opt-in before `posix_spawn`; its oracle inspects the constructed
environment while the parent value is set.  An explicit requested gate also rejects an
internal process worker and any team other than exactly 192 outer workers with
one inner thread.  Independently, a successful integrated batch arms an
implicit opportunity even when standalone permission is false and the opt-in is
absent.  That implicit path enables affinity only for an auto-selected 192x1
optimizer stage; smaller or otherwise different later optimizer topologies fall
back to normal CREST scheduling without error.  The readiness marker persists
across those stages and the next process-MTD batch entry/failure state clears it.

`optimizer_affinity_fortran_driver.F90` reproduces the production OpenMP
barrier/single/task control flow.  It requires zero tasks after either an
injected bind failure or a team-size mismatch, requires all tasks before an
injected post-restore error, and then requires a clean second region.  Both
drivers use actual per-thread Linux affinity and require full-mask restoration.

`check_optimizer_affinity_callers.py` scans the complete Fortran source tree,
requires the exact inventory of 14 direct `crest_oloop()` calls and five
transitive `crest_refine()` propagation calls, and rejects any site whose next
executable statement is not a typed-status guard containing cleanup-only work
and an immediate return. This prevents caller task, output, sorting, renaming,
or other postprocessing after a returned affinity failure.

Run both regular and AddressSanitizer/UndefinedBehaviorSanitizer variants into
a new evidence directory:

```sh
bash prototype_tests/run_optimizer_affinity_oracles.sh /absolute/new/output
```

The runner records separate stdout and stderr for every check, compiler, and
regular/sanitized executable. It hashes every local source file scanned or
compiled by the oracles and writes a complete SHA-256 manifest of the evidence
directory. The combined V9 build harness must generate a fresh full-tree
manifest with `compute_source_tree_sha256_v1.sh`; the V7-specific manifest
generator is deliberately not carried forward.

`nci_wall_lifecycle_driver.F90` proves that only the marked automatic wall is
detached, user constraints retain their order and data, the detached wall loses
its shallow process-local freeze pointer, and restoration rebinds it to the
current calculator freeze mask.

`fragment_string_copy_oracle.F90` and
`fragment_string_real_copy_oracle.F90` exercise the boxed scalar string
representation, real `calcdata%add` paths, nested copies, heap churn, and deep
copy independence. The capsule driver is the boxed-fragment version and must
still reproduce the exact capsule-v5 bytes.

These small-host oracles bypass only the production 192-CPU/Slurm discovery
precondition.  They do not replace the required full-node 48x4 child-mask,
192-worker optimizer-mask, science, or performance gates.

## Process-MTD worker-affinity oracle

`mtd_worker_affinity_oracle.c` is a two-CPU `posix_spawn`/`exec` regression for
the process-worker affinity boundary.  The driver selects two CPUs from its
real allowed mask without calling OpenMP, then execs a child with two explicit
singleton `OMP_PLACES` and `OMP_PROC_BIND=close`.  The child invokes the
production `crest_mtd_prepare_worker()` routine and immediately requires the
primary Linux-thread mask to be exactly the first singleton place.  The old
implementation deterministically fails here because it leaves the primary mask
as the full two-CPU worker set.

The fixed preparation applies only `{cpu0}` to the primary.  It never applies
the full worker set; the immediately following OpenMP team is the point-of-use
proof for both CPUs.  The child then forms one two-thread team and independently requires, for
each thread, agreement among the OpenMP singleton place CPU,
`sched_getaffinity()`, and `sched_getcpu()`, plus unique places and CPUs.  It
uses absolute Linux CPU IDs selected from the effective cpuset and needs no
CREST chemistry, large allocation, or Slurm submission.  An environment with
fewer than two allowed CPUs is reported explicitly with skip status 77.

This C oracle validates the production C affinity probe against independent
syscalls, but it does not directly call the private Fortran
`verify_worker_openmp_binding()` routine.  A clean C/Fortran link and the later
small integrated process-MTD run provide that separate coverage.

The runner also compiles the same independent oracle against the exact
optimizer-fallback-only source tarball.  That negative control is accepted only
when it fails at the widened-primary-mask assertion; any unrelated failure is
rejected.  After implementation review and explicit authorization to compile
and execute the oracle, run it from the source root with a new output directory:

```sh
bash prototype_tests/run_mtd_worker_affinity_oracle.sh \
  /absolute/new/output \
  /absolute/crest-process-mtd-dynamic-l3-optfallback-20260810.tar.gz
```
