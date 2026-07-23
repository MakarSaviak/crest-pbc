# CREST 3.1 host–guest workflow: build environment and external-rerank restart implementation

## 1. Purpose

This document is a reproducible implementation handoff for the CREST 3.1 host–guest workflow used for finite MOF–chromophore sampling. It covers:

1. how to restore and compile the exact CREST 3.1 source in the present constrained environment;
2. which already completed code should be used as the baseline;
3. how to implement the staged external-MLIP reranking restart;
4. how TOML-defined GFN-FF fragments must participate in the restart topology check;
5. how to test the implementation without assuming a frozen host.

The external MLIP calculation is a **ranking operation**, not a CREST energy or geometry replacement. After the first CREST iteration, an external workflow maps each finite-cluster guest pose into the periodic MOF, relaxes the guest with an MLIP, and ranks final MLIP energies. The selected restart structure passed back to CREST is the **original finite GFN-FF structure** corresponding to MLIP rank 1, not the periodic MLIP-relaxed structure.

---

## 2. Authoritative source baseline

### 2.1 Downloadable snapshot

Use the complete bundle:

```text
crest-3.1-host-guest-workflow-69d4fa9.tar.gz
```

Expected checksum:

```text
720e7941ddfb0fe92e289f0b0b9c7a79dc995eef7ddf811f5e6370a4168a9d46
```

The archive contains:

```text
source/                                      complete buildable source tree
crest-host-guest-workflow-31.git.bundle      Git history and branch
SNAPSHOT_MANIFEST.md                         provenance
```

The `source/` directory is the safest build tree because it contains the populated nested subprojects required by CMake. The Git bundle is the authoritative history/patch source. A Git clone reconstructed only from the bundle may lack nested subproject files that were not tracked by the top-level repository; in that case, use the complete `source/` tree for compilation and the Git clone for commits/diffs.

### 2.2 GitHub publication

Repository:

```text
MakarSaviak/AUGUR-MLIP
```

Completed Phase 1–5 branch:

```text
feature/host-guest-workflow-31
```

Draft pull request:

```text
PR #6 — Publish CREST 3.1 host–guest workflow phases 1–5
```

Published local source head represented by the patch series:

```text
69d4fa9fb522916037af32fc6d9875ea03fd8090
```

Baseline upstream CREST snapshot:

```text
pprcht/crest experimental
904c8b6a26b07b66705eeba490a2989e0c1f6e11
```

Completed local commits:

```text
838b985  separate MD/MTD thread budget (-TMD)
5c175dd  ordinary COM MTD bias and current-reference reseeding
beede30  multi-input NCI placements
b01b35a  reduced frozen-host ANC optimizer
69d4fa9  process-isolated ensemble scheduler
```

Do not start from the unfinished historical external-MLIP importer. The restart task described below replaces that more complicated design.

---

## 3. Current local environment

The working environment used for compilation and tests has approximately:

```text
usable CPU quota: 4 logical CPU equivalents
memory limit:     4 GiB
swap:             none
GPU:              none
C compiler:       GNU GCC 14.2
C++ compiler:     GNU G++ 14.2
Fortran compiler: GNU GFortran 14.2
build system:     CMake + Ninja
BLAS/LAPACK:      system libraries
OpenMP:           enabled
```

Because there is no swap, avoid unrestricted parallel compilation. `ninja -j2` is the safe default. `-j3` may work, but `-j2` leaves enough memory for large Fortran module compilations and linking.

---

## 4. Restore the source

### 4.1 Verify and extract the archive

```bash
sha256sum -c crest-3.1-host-guest-workflow-69d4fa9.tar.gz.sha256

tar -xzf crest-3.1-host-guest-workflow-69d4fa9.tar.gz
cd crest-3.1-host-guest-workflow-69d4fa9
```

### 4.2 Restore the Git history

```bash
git clone crest-host-guest-workflow-31.git.bundle crest-git
cd crest-git
git switch feature/host-guest-workflow-31
git rev-parse HEAD
```

Expected head:

```text
69d4fa9fb522916037af32fc6d9875ea03fd8090
```

Create a separate implementation branch:

```bash
git switch -c feature/external-rerank-restart-31
```

### 4.3 Recommended two-tree workflow

Use:

```text
crest-git/     Git authority: edits, diffs, tests added to version control
source/        complete build authority: populated nested subprojects
```

After editing files in `crest-git/`, copy only changed top-level CREST files into the corresponding locations under `source/` before compiling. Do not replace entire subproject directories from the Git clone.

---

## 5. Configure and build CREST 3.1

From the complete source tree:

```bash
cmake -S source -B build-external-rerank -G Ninja
ninja -C build-external-rerank -j2
```

Capture a resource report when benchmarking the build:

```bash
cd build-external-rerank
/usr/bin/time -v ninja -j2 \
  > ../build-external-rerank.log \
  2> ../build-external-rerank.time
```

Check the executable:

```bash
./build-external-rerank/crest --version
```

Check whether an incremental rebuild is needed:

```bash
ninja -C build-external-rerank -n
```

Expected result after a complete build:

```text
ninja: no work to do.
```

### 5.1 Build-system failure recovery

If Ninja was forcibly interrupted and reports an internal dyndep assertion, do not repeatedly restart the damaged build directory. Create a new build directory:

```bash
rm -rf build-external-rerank-clean
cmake -S source -B build-external-rerank-clean -G Ninja
ninja -C build-external-rerank-clean -j2
```

### 5.2 Unit tests

List tests:

```bash
ctest --test-dir build-external-rerank -N
```

Run the focused restart topology suite:

```bash
./build-external-rerank/test/crest-tester external_restart
```

Run earlier host–guest regressions:

```bash
test/host_guest_workflow/phase1_tmd_threads/run_phase1_checks.sh \
  ./build-external-rerank/crest

test/host_guest_workflow/phase2_com_reseed/run_phase2_checks.sh \
  ./build-external-rerank/crest

test/host_guest_workflow/phase3_multi_input_nci/run_phase3_checks.sh \
  ./build-external-rerank/crest

test/host_guest_workflow/phase6_external_rerank_restart/run_phase6_checks.sh \
  ./build-external-rerank/crest

tools/host_guest_workflow/tests/run_scheduler_smoke.sh
```

---

## 6. Target user workflow

### 6.1 First invocation

The scientific calculation should be described in TOML. The TOML contains the original input structure, calculation method, dynamics settings, GFN-FF fragments, atom selections, and target MTD iteration count.

Conceptual example:

```toml
runtype = "nci_search"
input = "initial-finite-host-guest.xyz"
external_rerank = true
mtd_iterations = 2
threads = 192
threads_md = 6

[[calculation.level]]
method = "gfnff"
fragments = ["1-808", "809-938"]
```

Run:

```bash
crest settings.toml
```

CREST performs:

```text
input parsing
optional initial optimization
trial MTD
MTD iteration 1
iteration-1 ensemble optimization
iteration-1 CREGEN
write iteration-1 archive
write staged restart checkpoint
exit normally before iteration 2
```

The checkpoint stage is:

```text
awaiting_external_rerank
```

### 6.2 External MLIP stage

Outside CREST:

1. Read the iteration-1 finite GFN-FF conformer ensemble.
2. Preserve a stable source-frame ID for every conformer.
3. Map only the guest from each finite cluster into the periodic MOF.
4. Relax the guest with the selected MLIP while using the intended host mobility constraints.
5. Rank by final periodic MLIP energy.
6. Identify the source CREST frame corresponding to MLIP rank 1.
7. Copy the **original finite GFN-FF full-system geometry** of that source frame to:

```text
crest-best-external.xyz
```

Do not write MLIP periodic energies into CREST’s GFN-FF ensemble. Do not use the periodic MLIP-relaxed structure directly as the CREST seed in the initial implementation.

### 6.3 Special restart invocation

Run:

```bash
crest crest-best-external.xyz --restart
```

No `--external-seed` option is required. The restart stage determines that the positional structure is the externally selected seed.

CREST must:

1. read `crest.restart`;
2. detect `stage = awaiting_external_rerank`;
3. load the original TOML path recorded by the checkpoint;
4. reparse the original scientific settings, including `fragments`;
5. override only the input geometry with `crest-best-external.xyz`;
6. skip ordinary preoptimization and trial MTD before validation;
7. validate the seed using fragment-aware topology;
8. restore iteration counters and iteration-1 archive state;
9. skip completed iteration 1;
10. run iteration 2 only;
11. collect iteration-1 and iteration-2 archives;
12. run final CREGEN and duplicate removal normally.

---

## 7. Required restart semantics

### 7.1 What the checkpoint stores

The special checkpoint should store workflow state, not duplicate every scientific setting:

```text
runtype
stage
completed MTD iteration
requested total MTD iterations
current nmetadyn
elowest
eprivious
iteration-1 archive path
original TOML settings path
```

The checkpoint does **not** need to store a second copy of the atom-by-atom fragment map. The authoritative fragment definition remains in the original TOML and is reparsed on restart.

### 7.2 Why the TOML path is needed

A normal new CREST process reconstructs `env` from defaults and the current command line. The historical lightweight `crest.restart` only stored stage counters and energies; by itself, it did not reconstruct calculation objects or TOML-defined fragments.

Therefore, the special checkpoint records the original TOML path. During:

```bash
crest crest-best-external.xyz --restart
```

CREST first reparses the saved TOML, then replaces its input geometry with the positional XYZ. This restores:

```text
calculation levels
GFN-FF fragments
charge and multiplicity
constraints and frozen atoms
RMSD/COM selections
NCI and wall settings
CREGEN settings
MTD settings
```

Runtime resources such as `-T` and `-TMD` may still be overridden explicitly on the restart command if desired.

### 7.3 Explicit restart requirement

An `awaiting_external_rerank` checkpoint must not be consumed automatically by launching CREST without `--restart`. This stage requires the explicit command:

```bash
crest crest-best-external.xyz --restart
```

This prevents the original input structure from being mistaken for the external seed.

---

## 8. Geometry and topology validation

### 8.1 Required checks

The externally selected structure must satisfy:

```text
same atom count
same element identities
same atom ordering
same finite/nonperiodic system type
same fragment-aware covalent topology as at least one iteration-1 frame
```

### 8.2 Checks that must not be imposed

Do **not** require:

```text
identical host coordinates
identical guest coordinates
identical full-system Cartesian coordinates
frozen host assumptions
```

The same restart mechanism must support a mobile host. If the original calculation froze the host, the reparsed constraints will freeze it during iteration 2. The seed validator itself must remain topology-based.

### 8.3 TOML fragments must participate in topology validation

For a TOML fragment definition such as:

```toml
fragments = ["1-808", "809-938"]
```

CREST/GFN-FF treats the listed groups as separate fragments and does not construct bonds between different fragment IDs. The topology validator must use the same partition.

The correct implementation is **not**:

```text
compute one full-system topology
then delete inter-fragment edges
```

That can fail when a close host–guest contact changes neighbour selection before masking.

The correct implementation is:

```text
resolve TOML fragment IDs with GFN-FF semantics
for each fragment independently:
    extract that fragment’s atoms
    compute its covalent topology
    compare reference and seed topology
ignore all inter-fragment contacts by construction
```

The implicit fragment `0` used for atoms omitted from explicit TOML groups must also be processed as its own fragment. Overlapping fragment definitions use the same last-definition-wins behavior as `gfnff_set_fragments()`.

If no user fragments are defined, compare the ordinary full-system topology.

### 8.4 Reference ensemble

The seed topology should be compared against the completed iteration-1 archive identified by `last_file` in the checkpoint. Accept the seed if it matches at least one iteration-1 frame. This avoids assuming that the original pre-sampling geometry is the only valid reference topology.

---

## 9. Exact implementation locations

### 9.1 `src/classes.f90`

Add workflow state to `systemdata`:

```text
input_settings_file
restart_requested
external_rerank
```

Update the deep-copy routine so queued/copied environments retain these fields.

### 9.2 `src/parsing/parse_maindata.f90`

Add TOML keys:

```text
external_rerank = true|false
mtd_iterations = integer
```

Aliases may include:

```text
mrest
maxrestart
iterations
```

Reject nonpositive iteration counts.

### 9.3 `src/confparse.f90`

Add CLI recognition for:

```text
--restart
--external-rerank
```

Before ordinary TOML discovery:

1. scan for explicit `--restart`;
2. read `crest.restart` if present;
3. if the stage is `awaiting_external_rerank` or `external_seed_loaded`, parse the saved TOML file;
4. set the positional non-TOML file as `env%inputcoords`;
5. disable preoptimization until the supplied geometry passes validation;
6. restore the target MTD iteration count from the checkpoint.

For an ordinary run, record the discovered TOML path in `env%input_settings_file`.

### 9.4 `src/restartlog.f90`

Extend `restart_data` with:

```text
target_mtd_iter
settings_file
```

Keep `write_restart_log()` backward-compatible by making the new arguments optional. Unknown keys must remain ignored for forward compatibility.

Add stage labels:

```text
awaiting_external_rerank
external_seed_loaded
```

### 9.5 `src/external_rerank_restart.f90`

Add a dedicated module containing:

```text
validate_external_rerank_seed()
fragment_aware_topology_equal()
resolve_gff_fragment_ids()
```

The module should:

- retrieve the TOML-reconstructed `gff_fragments` from calculation levels;
- require consistent fragment definitions across levels that define them;
- reproduce GFN-FF fragment parsing semantics;
- compare each fragment topology independently;
- compare the external seed against iteration-1 archive frames;
- never compare host coordinates for equality.

Add the module to:

```text
src/CMakeLists.txt
src/meson.build
```

### 9.6 `src/algos/search_conformers.f90`

At restart detection:

- require explicit `--restart` for `awaiting_external_rerank`;
- distinguish the special external restart from ordinary restart stages.

Before trial MTD:

- validate the supplied seed;
- write `external_seed_loaded` immediately after validation;
- skip trial MTD for all restarts.

Inside the MTD loop:

- skip iterations `<= completed_iteration` for the external stages;
- seed the next MTD from the new `env%ref` supplied by the positional XYZ.

After completing iteration 1 in an external-rerank run:

- finish iteration-1 cleanup/state updates exactly as ordinary CREST does;
- write `stage = awaiting_external_rerank`;
- record the iteration-1 archive, target iteration count, and TOML path;
- return normally before iteration 2.

### 9.7 `src/cleanup.f90`

CREST normally removes `.cre_*` files during the global process cleanup. That behavior is correct after a completed ordinary run, but it would destroy the iteration-1 archive immediately after the intentional external-rerank pause. In `custom_cleanup(env)`, preserve `.cre_*` while `env%external_rerank` is active. Final `collectcre` remains responsible for consuming/removing the archives after all requested iterations are complete.

This is not optional: the checkpoint can only be resumed and the final cross-iteration CREGEN can only be correct if the iteration-1 archive survives process 1.

### 9.8 Tests

Add:

```text
test/test_external_restart.F90
```

Register suite:

```text
external_restart
```

in:

```text
test/CMakeLists.txt
test/meson.build
test/main.f90
```

Required unit cases:

1. two fragments preserve their internal topology while moving close enough to create inter-fragment contacts: fragment-aware comparison passes;
2. the same structures without fragments: full topology comparison detects the difference;
3. an internal bond change inside one fragment: fragment-aware comparison fails;
4. fragment range strings resolve to expected per-atom IDs.

Add an end-to-end shell test that:

1. runs a very short two-iteration TOML calculation with `external_rerank = true`;
2. confirms the first process stops at `awaiting_external_rerank` after iteration 1;
3. creates an external seed from the iteration-1 archive;
4. runs `crest external.xyz --restart`;
5. confirms iteration 1 is skipped and iteration 2 runs;
6. confirms final collection reaches `stage done` after using both iteration archives;
7. rejects a seed with an intrafragment topology change before iteration 2 begins.

The mobile-host rule is exercised directly in the Fortran unit suite: an entire fragment is translated while its internal topology is preserved, and the fragment-aware comparison must pass. This isolates the restart validator from unrelated downstream sampling behavior of a tiny smoke-test molecule.

---

## 10. State machine

```text
new run
  |
  v
iteration_1_running
  |
  v
iteration_1_complete
  |
  v
awaiting_external_rerank
  |
  |  crest crest-best-external.xyz --restart
  v
external_seed_validated
  |
  v
external_seed_loaded
  |
  v
iteration_2_running
  |
  v
all_requested_iterations_complete
  |
  v
collect iteration archives
  |
  v
final CREGEN
  |
  v
done
```

If the process is interrupted after `external_seed_loaded`, rerunning the same restart command must not launch iteration 1 again. The checkpoint’s completed-iteration count remains authoritative.

---

## 11. Acceptance criteria

The implementation is acceptable when all of the following hold:

- The Phase 1–5 tests remain green.
- `crest settings.toml` pauses normally after iteration 1 when `external_rerank = true`.
- `crest crest-best-external.xyz --restart` automatically reparses the original TOML.
- TOML GFN-FF fragments are present before GFN-FF initialization in iteration 2.
- The topology check ignores inter-fragment contacts by computing fragment topologies independently.
- A mobile host is accepted when its intrafragment topology is unchanged.
- An intrafragment bond change is rejected before MTD or optimization.
- Iteration 1 is not rerun.
- The target number of iterations is not increased by restart.
- Iteration-1 archives remain available for final CREGEN.
- CREST and MLIP energy scales remain separate.
- No MLIP-relaxed periodic structure is silently substituted for a finite GFN-FF structure.

---

## 12. Git workflow after validation

From the Git-authority tree:

```bash
git status --short
git diff --check
git add \
  src/classes.f90 \
  src/cleanup.f90 \
  src/parsing/parse_maindata.f90 \
  src/confparse.f90 \
  src/restartlog.f90 \
  src/external_rerank_restart.f90 \
  src/algos/search_conformers.f90 \
  src/CMakeLists.txt \
  src/meson.build \
  test/CMakeLists.txt \
  test/meson.build \
  test/main.f90 \
  test/test_external_restart.F90 \
  test/host_guest_workflow/phase6_external_rerank_restart

git commit -m "Add external-rerank staged restart"
```

Push to a branch separate from the Phase 1–5 publication:

```text
feature/external-rerank-restart-31
```

The PR should target the Phase 1–5 branch while under development, or `main` only after the Phase 1–5 patch publication is incorporated into the chosen GitHub history.

Do not mix the abandoned full MLIP-ensemble importer into this commit.


---

## 13. Implementation completed in the present session

### 13.1 Local branch and source trees

The implementation was made on:

```text
feature/external-rerank-restart-31
```

with Phase 1–5 commit `69d4fa9fb522916037af32fc6d9875ea03fd8090` as its parent.

The Git-authority tree used for diffs and commits was:

```text
/mnt/data/crest_restart_doc_work/repo
```

The complete source tree used for compilation was:

```text
/mnt/data/crest_restart_doc_work/full
```

Only changed top-level CREST files were synchronized from the Git-authority tree into the complete build tree. This avoids losing populated nested subproject content while keeping the commit history clean.

### 13.2 Implemented behavior

The implementation now provides the following staged workflow:

```bash
# Process 1
crest input.toml

# External MLIP ranking writes/copies the selected original finite frame

# Process 2
crest crest-best-external.xyz --restart
```

The first process runs iteration 1 and exits normally at:

```text
stage awaiting_external_rerank
mtd_iter 1
target_mtd_iter 2
settings_file input.toml
```

The iteration-1 `.cre_0.xyz` archive is preserved across process cleanup.

The second process:

1. recognizes the explicit `--restart` and special checkpoint stage;
2. reparses the original TOML before calculator construction;
3. replaces only the positional geometry with `crest-best-external.xyz`;
4. obtains `fragments = [...]` from the reconstructed GFN-FF calculation object;
5. validates atom count, elements/order, and fragment-internal covalent topology against the iteration-1 archive;
6. imposes no coordinate-equality requirement on either host or guest;
7. skips iteration 1;
8. runs iteration 2;
9. collects the iteration archives and reaches `stage done`.

No fragment list is duplicated into the custom restart checkpoint. The checkpoint stores the original TOML path; the normal TOML parser remains the single source of truth for fragments and the other scientific settings.

### 13.3 Files changed

Implementation files:

```text
src/classes.f90
src/cleanup.f90
src/confparse.f90
src/parsing/parse_maindata.f90
src/restartlog.f90
src/external_rerank_restart.f90
src/algos/search_conformers.f90
src/CMakeLists.txt
src/meson.build
```

Test files:

```text
test/test_external_restart.F90
test/CMakeLists.txt
test/meson.build
test/main.f90
test/host_guest_workflow/phase6_external_rerank_restart/input.toml
test/host_guest_workflow/phase6_external_rerank_restart/struc.xyz
test/host_guest_workflow/phase6_external_rerank_restart/run_phase6_checks.sh
```

### 13.4 Build result

A clean CMake/Ninja build completed successfully:

```text
build directory: /mnt/data/crest_restart_doc_work/build-external
targets:         1598/1598
parallelism:     ninja -j2
wall time:       4 min 49.68 s
CPU utilization: 199%
peak RSS:        1,184,380 KiB
exit status:     0
```

The built executable reports CREST `3.1.0`.

### 13.5 Tests completed

The new Fortran suite passes all three cases:

```text
external restart fragment topology       PASSED
external restart internal bond change    PASSED
external restart fragment resolution     PASSED
```

The focused Phase 6 shell integration test passes:

```text
iteration 1 executes
awaiting_external_rerank checkpoint is written
.cre_0.xyz survives process cleanup
original TOML and two GFN-FF fragments are restored
coordinate equality is not required
iteration 1 is skipped on restart
iteration 2 executes
valid continuation reaches stage done
intrafragment topology corruption is rejected before iteration 2
```

The complete CMake/CTest matrix passed: tests 1–36 passed in the first run, and tests 37–68 passed in the continuation run. After the final source cleanup, all 15 CREST-specific suites (`crest/tblite` through `crest/external_restart`) were rerun and passed again.

The host–guest workflow regressions also pass with the final executable:

```text
phase1 -TMD checks passed
phase2 COM/reseed checks passed
phase3 multi-input NCI checks passed
phase6 external-rerank restart checks passed
process scheduler smoke test passed
```

### 13.6 Tested fragment semantics

For a TOML calculation containing:

```toml
[[calculation.level]]
method = "gfnff"
fragments = ["1-3", "4-6"]
```

restart validation constructs topology separately for fragment 1 and fragment 2. Inter-fragment proximity does not create a topology mismatch. Translating a complete fragment is therefore accepted by the validator, while stretching an O–H bond inside one fragment is rejected.

This is the intended mobile-host behavior: coordinates may change; fragment membership and intrafragment covalent connectivity must not.
