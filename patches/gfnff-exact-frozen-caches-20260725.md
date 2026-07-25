# Exact frozen-host GFN-FF caches

This incremental patch applies after `patches/gfnff-active-force-loops-20260725.patch`.

## Scope

For repeated GFN-FF energy/gradient calls with a fixed host and a smaller active region, the patch adds persistent work arrays and exact caches for quantities that are invariant while the frozen coordinates remain unchanged:

- frozen–frozen distances;
- frozen-prefix distance and D3-list setup;
- raw frozen–frozen coordination-number contributions;
- fixed-form frozen-only nonbonded and bonded repulsion energies;
- all-frozen angle and torsion energies;
- dynamic lists for the remaining terms.

The patch does not freeze, approximate, or reduce the update frequency of EEQ charges, active–frozen or active–active CN terms, D3 coefficients or energies, CN-dependent SRB bonds, HB/XB terms, or any active force.

Caches are invalidated on changes to the freeze mask, frozen coordinates, topology dimensions, or relevant thresholds. Non-prefix masks use a general exact fallback.

## Validation

The benchmark system contained 938 atoms, with atoms 1–808 frozen and atoms 809–938 active. No structures, trajectories, chemical identities, or project-specific names are published here.

Across 12 separated structures, cache versus active-force parent maxima were:

- total/decomposed energy difference: `1.71e-13 Eh`;
- EEQ charge difference: `2.89e-15 e`;
- active-gradient component difference: `1.11e-16 Eh/bohr`.

The small energy difference is due to floating-point summation grouping of static scalar terms. The active gradient agrees to numerical precision. Comparisons at 1, 2, and 4 OpenMP threads stayed within the parent implementation's own cross-thread reduction variation.

A deterministic constrained relaxation gave all three compared builds the same 23 force calls, final energy, final gradient norm, and byte-identical optimized coordinates.

Seven 200-step MTD runs completed normally. All retained the fixed host and produced stable, overlapping energy and structural distributions. The stochastic trajectories were not expected to match coordinate-by-coordinate.

## Performance

Sandbox, one CPU core:

| Test | Cache | Active-force parent | Previous exact-frozen | Cache gain vs active parent |
|---|---:|---:|---:|---:|
| Repeated force calls | 98.76 ms/call | 104.73 ms/call | 110.06 ms/call | 5.70% |
| Deterministic relaxation | 7.55 s | 7.65 s | 8.09 s | 1.31% |
| 200-step MTD mean | 24.50 s | 24.89 s | 28.21 s¹ | 1.55% |

¹ One previous-exact run; cache and active-force values are three-run means.

HPC-node benchmarking is still required because CPU-frequency behavior and OpenMP scaling are hardware-dependent.

## Tests

- complete configured CTest matrix: 9/9 passed;
- GFN-FF test subset passed at 1, 2, and 4 OpenMP threads;
- frozen-host regression covers cache creation, reuse, frozen-coordinate invalidation, mask invalidation, non-prefix masks, and freeze-mask removal.
