# Exact frozen-host MTD acceleration validation

## Scope

This follow-up patch accelerates GFN-FF energy/gradient calls used by CREST when atoms are frozen.

CREST passes its atom-freeze mask into the GFN-FF calculator. GFN-FF then omits only direct force accumulation on frozen endpoints for:

- nonbonded repulsion;
- bonded repulsion;
- direct electrostatic pair gradients;
- direct D3 pair gradients.

The calculation still evaluates the complete potential energy and all quantities that can influence active atoms, including:

- coordination numbers and their derivatives;
- the full EEQ charge solution;
- electrostatic energy and CN-dependent electrostatic gradients;
- D3 energy, coefficients, and CN-dependent D3 gradients;
- SRB bonds, angles, torsions, bonded ATM, HB/XB, and other GFN-FF terms;
- all active-active and active-frozen forces.

The optimized path therefore changes only frozen reaction-force components that CREST would immediately zero before dynamics. No force-field parameter, cutoff, convergence criterion, update frequency, charge, CN, energy, or active force is approximated.

## Source provenance

- Original validated host-guest source head: `092377fe68dd6cf493c98c377aca8861213d859c`
- Exact frozen-gradient implementation: `74c26b4`
- Full validated source head: `932f6b0b810b1d3a62204bb181552a9de7ef27f4`

## Frame-by-frame equivalence

Validation used the real 12-frame, 938-atom AM03@MOF trajectory included with the source bundle. The production mask freezes atoms 1-808 and leaves AM03 atoms 809-938 active.

At one OpenMP thread, for all 12 frames:

- total energies were bitwise identical;
- EEQ charges were bitwise identical;
- every decomposed GFN-FF energy component was bitwise identical;
- every active-atom gradient component was bitwise identical;
- the complete gradient after CREST's normal freeze operation was bitwise identical.

Additional masks were tested successfully:

- allocated all-active mask;
- sparse frozen-host mask;
- host plus selected guest atoms frozen;
- nearly all atoms frozen.

At two and four OpenMP threads, the maximum differences were no larger than the ordinary OpenMP reduction variation observed when repeating the unoptimized calculation:

- maximum energy difference: `2.8421709430404007e-14 Eh`;
- maximum active-gradient difference: `1.3877787807814457e-17 Eh/bohr`;
- maximum decomposed-component difference: `2.8421709430404007e-14 Eh`.

## Performance

### Repeated 938-atom GFN-FF force loop

Representative paired one-thread measurements:

- baseline: `17.08 s`;
- optimized: `15.14 s`;
- improvement: `11.4%`.

A longer development measurement gave:

- baseline: `30.28 s`;
- optimized: `26.30 s`;
- improvement: `13.1%`.

### Standalone CREST MTD benchmark

Five paired one-thread runs used the real 938-atom AM03@MOF system, 808 frozen host atoms, explicit host/guest fragments, RMSD metadynamics, and 20 energy/gradient steps.

Mean wall times:

- baseline: `7.004 s`;
- optimized: `6.624 s`;
- improvement from means: `5.425%`;
- mean paired improvement: `5.337%`.

The smaller complete-MTD gain is expected because CREST performs setup and dynamics work outside the GFN-FF force routine.

## Regression tests

- complete configured CTest matrix: `9/9 passed`;
- GFN-FF suite: `6/6 passed` at 1, 2, and 4 OpenMP threads;
- new regression: `GFN-FF exact frozen gradient`.

The regression verifies equal energy, equal active gradients, zero returned frozen gradients, and equality after CREST applies its normal freeze operation.

## Deferred experiments

Persistent scratch allocation, frozen-pair distance caching, and broader static-term caching showed approximately `3.5-4.1%` combined potential in development testing. They are deliberately not enabled in this production patch because the simpler exact-gradient change is better isolated and already clears the end-to-end acceptance threshold.

The complete downloadable bundle includes the deferred experimental patch separately for later work.
