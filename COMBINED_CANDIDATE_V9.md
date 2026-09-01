# HISTORICAL — DOES NOT DESCRIBE THE CURRENT SOURCE

This archived V9 provenance document predates the current generic process-MTD
implementation. Its V8/V9 identity claims, fixed-system examples, capsule-V4
statements, and opt-in wording are historical only. Do not use it as a release
or validation guide; see `CURRENT_PROCESS_MTD_STATUS.md` instead.

# Combined fragment-box and optimizer-affinity candidate V9

Status: documentation-corrected source-composition candidate only. It is not
a production build and has not passed the final-executable scientific,
full-node, or performance gates.

V9 is a documentation-only successor to frozen V8 tree SHA-256
`b007fb57354cb5d47b4db9a53e95bad418f669e0d7d19bad64dc4fbe35c1c226`.
All production and prototype-test sources are byte-identical to V8. The stale
inherited claim that this source was already final, and the invalid claim of
an independently enforced single BLAS thread per caller, were removed before
any combined clean build.

## Frozen inputs

- Optimizer-affinity V7 full-tree SHA-256:
  `ec2714f4c5054f46801b5bcdf03b3e89a21c9e28c1bba866fd4939a5c487fa19`.
- Boxed-fragment manifest SHA-256:
  `4e0a6b7521abc65847107c3a515d45c1cf861dd37ee993679a32f29bf58aae91`.
- Common frozen NCI-walls V5 predecessor: the preserved source under
  `builds/C_hbreset_release_v6_process_isolated_mtd_scheduler_always_real64_eeq_fitde_historysafe_statusprop_nciwalls_v5/source`.

## Exact composition

The production composition originated from V7 without its generated
`PROTOTYPE_EVIDENCE` directory. Five fragment-only production files are byte-identical to the
boxed-fragment branch. `src/mtd_process_scheduler.F90` is a conflict-free
three-way merge of V7 (ours), V5 (base), and boxed-fragment V1 (theirs), with
SHA-256
`91bfe1ba35c37fe55d0750c428115f31ef7a57eacb96e2dc3a884f87df653d6d`.
The combined candidate intentionally excludes experimental D3 changes.

Relative to V5, exactly 14 production files may differ: eight affinity-only
files, five fragment-only files, and the composed scheduler. The expected
production diffstat is 622 insertions and 32 deletions. The complete GFN-FF
tree and the 60-line optimizer scientific task body must remain byte-identical
to V5. The capsule-v3 wire fixture must remain byte-identical with SHA-256
`6df219ca2802b5a4dff73019e7d3f28635be29eb5fb8bde60980ff3aeef4ab40`.

## Required gates

1. Verify the exact source decomposition and create a new full-tree manifest.
2. Rerun regular and ASan/UBSan affinity, fragment-copy, real-copy,
   capsule-v3, process-option, and caller-hardening oracles.
3. Build from empty source/build directories with GCC/GFortran 12.3,
   x86-64-v3 flags, OpenMP, and pinned NUM_THREADS=256 OpenBLAS; preserve all
   compiler, linker, binary, and dependency evidence.
4. Run parser, fragment-order, capsule, external-rerank, Fit-De/history, and
   T1/T2/T4 final-executable scientific gates.
5. Prove full-node stage-specific affinity and mask restoration, then run the
   integrated 8x6 process-MTD gate and audit every child and downstream stage.
6. Run the fixed 54,300-frame screen and only then a one-iteration
   production-like gate. Stochastic MTD timings are compared as balanced
   distributions, never by trajectory identity.
