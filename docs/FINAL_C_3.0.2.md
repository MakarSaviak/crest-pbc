# CREST 3.0.2 final C candidate

This branch contains the final C implementation developed for a frozen-host
host–guest GFN-FF workflow. It retains the exact frozen-host EEQ decomposition,
the multi-input NCI workflow, COM metadynamics bias, custom reseeding, and the
MTD synchronization cleanup.

The final hot-path fixes:

- copy a frozen mask only when it is allocated, resized, or changed;
- parse explicit fragment assignments only when GFN-FF topology is initialized
  or refreshed;
- validate frozen coordinates with a scaled floating-point tolerance without
  allocating full-size temporary arrays; and
- remove the unconditional ANC optimizer diagnostic-matrix dump.

## Validation summary

The production-shaped validation system had 938 atoms: an 808-atom frozen host
and a 130-atom active guest. No molecular coordinates, trajectories, or private
calculation outputs are included in this repository.

Against the cache-only full-EEQ control built with identical compiler and
linkage settings, final C agreed to double-precision roundoff across five input
frames: maximum energy-component difference `2.27e-13 Eh`, maximum live atomic
charge difference `8.77e-15 e`, and maximum active-gradient-component
difference `1.73e-16 Eh/Bohr`.

Against the official HPC CREST 3.0.2/GFN-FF implementation running strictly
one software thread, final C agreed within `2.84e-13 Eh` in total energy,
`4.08e-16 Eh/Bohr` in active-gradient components, and `1.19e-14 e` in live
atomic charges across the same five fixed-coordinate inputs. Both builds
converged all five structures under crude and tight optimization, with the
frozen host unchanged exactly.

Two reversed-order, full production-shaped comparisons completed normally.
Final C's MTD stage took 76:41 and 77:41, compared with 83:25 and 84:22 for the
control: a 7.99% lower mean MTD wall time. Total C runtimes were 2:24:21 and
2:19:53; control runtimes were 2:30:01 and 3:10:11. The total-runtime spread
was dominated by a stochastic optimization long tail, so the two-repetition
mean is not presented as a portable effect size.

High-concurrency builds must use a BLAS implementation configured for at least
the number of concurrent callers. The validated 192-worker binaries used a
private OpenBLAS build with `MAX_THREADS=256` and one BLAS thread per caller.

## Reference scope

Runtime comparisons concern stochastic workflow timing, not trajectory or
endpoint identity. Correctness comparisons use fixed coordinates and report
energies, charges, and gradients separately.
