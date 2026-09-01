# Validation record — 2026-08-26

## Inputs and build

- Baseline source archive SHA-256:
  `20bdc4f466bdd4b02891aabd47f93d1b4e99fbfeb2c8f6e4f68ff2d2764acbd4`
- Supplied OpenBLAS archive SHA-256:
  `0c9a58a891c19d8607783ac5e35a069aed91ba89cf0575332d5f5f6092f7eeca`
- Toolchain: GCC/GFortran 14.3.0.
- OpenBLAS runtime: `0.3.24`, `USE_OPENMP`, `MAX_THREADS=256`.

The full feature build used installed tblite, toml-f, mctc-lib, s-dftd3,
dftd4, multicharge, and test-drive from the existing Conda prefix. GFN-FF,
GFN0, and lwONIOM were built from the vendored source subprojects because the
installed GFN-FF/GFN0 package configurations do not expose compiler-compatible
Fortran module files. The executable links the supplied OpenBLAS library.

For this OpenMP OpenBLAS build, `OPENBLAS_NUM_THREADS=1` alone does not change
the initial OpenBLAS OpenMP team count. CREST invokes `openblasset(1)` before
worker calculator work; both the C and Fortran OpenBLAS setters were checked
against this library and changed `openblas_get_num_threads()` from 192 to 1.

## Passing validation

- Clean full-feature CMake build with the supplied OpenBLAS-256 library.
- CTest: `pvol`, `gfnff`, `gfn0`, `gfn0occ`, `CN`, `optimization`, and the
  new `mtd-capsule-v5` regression.
- `mtd-capsule-v5`: AM03-shaped 938-atom regression plus 17- and 47-atom
  generic cases; V5 byte round trip; corruption gates; and explicit V4-header
  rejection.
- Dynamic L3 planner including `K=8` and child `BLAS=1` environment contract.
- Worker-affinity child oracle.
- Optimizer-affinity C, Fortran, sanitizer, and caller-propagation oracles.

## Known non-passing check and remaining release gates

The full CTest run has one failure: `crest/tblite` GFN2-xTB ALPB-water expects
`-23.983234299793384` but the installed tblite dependency stack returns
`-23.976134262859187`. Its other tblite cases pass, as do the GFN-FF process
dependencies. This numerical reference mismatch is not changed or hidden by
this source patch and must be resolved or accepted against a pinned dependency
stack before a full release.

No real small non-AM03 multi-trajectory GFN-FF NCI/RMSD-MTD fixture or AM03
production harness is included in this source tree. Those two end-to-end runs
remain mandatory release gates. This is a source candidate, not a promoted
production build.
