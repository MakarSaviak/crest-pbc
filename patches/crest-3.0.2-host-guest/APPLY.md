# CREST 3.0.2 host-guest workflow patch set

This directory publishes the validated CREST 3.0.2-v2 host-guest workflow and the exact frozen-host MTD acceleration developed for the AM03@MOF workflow.

## Baseline and validated heads

- CREST repository baseline: `cfdc301f759686b0fd66ced63b5ddbd6c693fa4f`
- Validated host-guest source head: `092377fe68dd6cf493c98c377aca8861213d859c`
- Validated exact-MTD source head: `932f6b0b810b1d3a62204bb181552a9de7ef27f4`
- GFN-FF remains the upstream CREST submodule; compatibility and acceleration changes are applied as separate verified patches.

## Apply

From a checkout of the baseline commit with this publication directory present:

```bash
git checkout cfdc301f759686b0fd66ced63b5ddbd6c693fa4f
git submodule update --init --recursive
bash patches/crest-3.0.2-host-guest/apply.sh
```

The script applies, in order:

1. the original split `crest-3.0.2-host-guest.patch` to the top-level CREST source;
2. `crest-exact-frozen-mtd-gradients.patch` to the top-level CREST source and tests;
3. `gfnff-host-guest-fragments.patch` to `subprojects/gfnff`;
4. `gfnff-exact-frozen-mtd-gradients.patch` to `subprojects/gfnff`.

Every patch is checked against its SHA-256 checksum before application.

## Host-guest workflow features

- separate `-TMD` dynamics/metadynamics thread budget;
- ordinary COM metadynamics bias and current-reference reseeding;
- multi-input NCI placement support;
- explicit host/guest fragmentation for the CREST 3.0.2 GFN-FF backend;
- staged external-rerank restart from the selected original CREST source frame.

## Exact frozen-host MTD acceleration

CREST passes the frozen-atom mask into GFN-FF. GFN-FF skips only direct reaction-force accumulation on frozen endpoints for selected pairwise terms. The full energy, EEQ solution, coordination-number machinery, dispersion model, CN-dependent gradients, and every force on active atoms remain unchanged.

The optimization does not alter:

- GFN-FF parameters or force-field equations;
- EEQ or CN update frequency;
- solver tolerances;
- cutoffs;
- active-atom forces;
- total or decomposed energies.

See `EXACT_FROZEN_MTD_VALIDATION.md` for the full equivalence and performance record.

## Validation summary

For 12 real 938-atom AM03@MOF frames with atoms 1-808 frozen, one-thread energies, charges, decomposed energies, active gradients, and post-freeze gradients were bitwise identical.

Observed performance:

- repeated GFN-FF force loop: approximately `11-13%` faster;
- five-run standalone CREST MTD mean: `5.425%` faster.

Regression results:

- complete configured CTest matrix: `9/9 passed`;
- GFN-FF suite: `6/6 passed` at 1, 2, and 4 OpenMP threads.

## Validated downloadable bundle

Bundle SHA-256:

```text
8533dd6b52edfd26fcc72e90268dbe122d286dcd53cf09e193a280855955228e
```

The bundle contains the complete source Git history, validated executable, build instructions, validation records, publication patches, and a separate non-applied experimental cache patch for possible later development.
