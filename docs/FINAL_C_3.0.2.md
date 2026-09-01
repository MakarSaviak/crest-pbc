# Historical release-v3 C evidence—not current V9 acceptance

This file records evidence for an earlier release-v3 C implementation of the
frozen-host host–guest workflow. It is retained as provenance only. The
combined fragment-box and optimizer-affinity V9 source is a new, source-only
candidate and is not final, production accepted, or covered by the historical
timings below.

The historical implementation retained exact frozen-host EEQ, multi-input
handling, COM metadynamics bias, custom reseeding, and synchronization
cleanup. Its hot-path changes copied the frozen mask only when needed, moved
fragment parsing to topology initialization/refresh, used allocation-free
scaled frozen-coordinate validation, and removed the unconditional ANCOPT
diagnostic-matrix dump.

## Historical validation summary

The private validation system had 938 atoms: an 808-atom frozen host and a
130-atom active guest. No molecular coordinates, trajectories, or private
calculation outputs are included in the source tree.

Against the cache-only full-EEQ control built with identical compiler and
linkage settings, the earlier release-v3 C build agreed to double-precision
roundoff across five input frames: maximum energy-component difference
`2.27e-13 Eh`, maximum live atomic-charge difference `8.77e-15 e`, and
maximum active-gradient-component difference `1.73e-16 Eh/Bohr`.

Against official HPC CREST 3.0.2/GFN-FF running one software thread, that
earlier build agreed within `2.84e-13 Eh` in total energy,
`4.08e-16 Eh/Bohr` in active gradients, and `1.19e-14 e` in live charges.
Both builds converged all five crude and tight optimizations while leaving the
frozen host unchanged.

Two reversed-order historical workflow comparisons completed normally. The
earlier C MTD stage took 76:41 and 77:41, versus 83:25 and 84:22 for its
control. Total runtimes were 2:24:21 and 2:19:53 for C, versus 2:30:01 and
3:10:11 for the control. These stochastic two-repetition timings are not a
portable effect size and do not validate V9.

## OpenBLAS correction

High-concurrency executables use the private OpenBLAS 0.3.24 build configured
with `MAX_THREADS=256`, preventing the earlier 128-slot metadata overflow.
That library is an OpenMP build: `OPENBLAS_NUM_THREADS=1` and the existing
setter are not independent one-thread kernel controls, because the stored
count follows the active OpenMP team. Therefore this source does not claim
“one BLAS thread per caller.” Final threading and affinity must be proved from
the exact V9 executable in allocated runtime gates.

## Scope

Stochastic workflow timing never implies trajectory or endpoint identity.
Correctness comparisons use fixed coordinates and report energy, charge, and
gradient differences separately. V9 still requires a clean build, scientific
oracles, full-node affinity proof, fixed 54,300-frame screen, and integrated
one-iteration validation before any production or final-candidate claim.
