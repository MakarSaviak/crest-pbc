# CREST 3.1 host-guest workflow extensions

This branch publishes the **implementation only** of the CREST 3.1
host-guest workflow.  It intentionally contains no production host-guest
geometry, trajectory, conformer ensemble, benchmark result, Slurm output, or
external-MLIP relaxation data.  The small XYZ files below `test/` are
synthetic regression fixtures only.

## Provenance

The branch starts from the CREST 3.1 experimental snapshot
`904c8b6a26b07b66705eeba490a2989e0c1f6e11` and contains the validated local
implementation series ending at:

```text
43d7b137910a69990db87c043f89ac49dcf0876e  Add external-rerank staged restart
```

It is deliberately separate from the 3.0.1-derived development branches in
this repository.  It neither rebases those branches nor imports unrelated
CREST 3.0.2 code.

## What this branch adds

- A separate `-TMD` thread budget for MD/MTD work.
- COM metadynamics with current-reference reseeding.
- Multi-input NCI starting placements.
- Reduced frozen-host ANC optimization, including direct active-coordinate
  Hessian construction and thin ANC-mode storage while retaining the ordinary
  non-frozen optimization path.
- A process-isolated ensemble scheduler, designed to run one independent
  single-thread CREST process per frame.
- A staged external-MLIP reranking restart: CREST pauses after a completed
  first MTD iteration; an external workflow ranks the finite-system frames;
  CREST resumes from the corresponding original finite GFN-FF geometry.

The external MLIP result is used only to choose a source frame.  Its energy
does not replace a CREST/GFN-FF energy, and a periodic MLIP-relaxed geometry
is not passed directly back to CREST.

## External-rerank restart workflow

1. Start a TOML-defined calculation with `external_rerank = true` and a
   positive `mtd_iterations` value.
2. CREST completes iteration 1, preserves its archive, writes
   `crest.restart` with stage `awaiting_external_rerank`, then exits normally.
3. Externally rank the iteration-1 frames while retaining each frame's stable
   source ID.  Write the selected frame's original finite GFN-FF coordinates
   as `crest-best-external.xyz`.
4. Resume explicitly:

   ```bash
   crest crest-best-external.xyz --restart
   ```

On restart, CREST reparses the original TOML settings, restores its fragment
definitions and constraints, validates the supplied seed against the
iteration-1 archive using fragment-aware topology, skips completed MTD
iterations, and performs the remaining work.  The validation permits a
mobile host; it does not require Cartesian identity with an earlier frame.

## Build and tests

Use a complete source tree with its populated nested subprojects:

```bash
cmake -S . -B build -G Ninja
ninja -C build -j2
ctest --test-dir build
```

The detailed build, restart, fragment-topology, and focused test instructions
are in [docs/host_guest_external_rerank_restart.md](docs/host_guest_external_rerank_restart.md).
The added end-to-end fixture is under
`test/host_guest_workflow/phase6_external_rerank_restart/`.

## Important operational boundary

This repository branch is code and synthetic tests.  Keep research inputs,
full trajectories, conformer structures, scheduler workspaces, raw logs, and
MLIP outputs in private storage.  Do not add them to commits or pull requests.
