# CREST 3.0.2 multi-input canonical/preoptimization fix

This public-safe patch fixes the multi-input NCI initialization failure that occurred when ordinary CREST preoptimization changed the first prepared structure before the multi-input loader validated the original positional ensemble.

## Required base

Apply the existing CREST 3.0.2 host–guest workflow first:

```bash
git checkout cfdc301f759686b0fd66ced63b5ddbd6c693fa4f
git submodule update --init --recursive
bash patches/crest-3.0.2-host-guest/apply.sh
```

Then apply this fix:

```bash
bash patches/crest-3.0.2-multiinput-canonical-preopt-fix/apply.sh
```

The patch checksum is:

```text
7433892e48df721f77b8978e915bc5336409253e5197cffafcf382c60d7caeac  multiinput-canonical-preopt-fix.patch
```

## Root cause

The ordinary CREST input path writes the transformed but unoptimized first structure to `coord`. During NCI preparation, `prepared_ref` may then be changed by geometry optimization. The previous multi-input implementation transformed raw frame 1 and compared it directly with this optimized `prepared_ref` using a strict Cartesian tolerance. The comparison therefore mixed a coordinate-frame transformation with real internal geometry relaxation and failed even for valid inputs.

Bypassing the check alone would also have been incorrect because it would rebuild frame 1 from the raw input and discard the optimized first seed.

## Correct behavior

The patch:

- treats `coord` as the authoritative transformed, unoptimized canonical reference;
- determines one proper rigid transformation from raw frame 1 to `coord` with CREST's quaternion least-squares RMSD routine;
- applies exactly the same rotation and translation to every supplied input frame;
- retains the optimized `prepared_ref` as production frame 1;
- preoptimizes only additional frames 2–N, matching the validated CREST 3.0.1 lifecycle;
- verifies a determinant near `+1` and retains strict RMSD and Cartesian postconditions;
- rejects non-finite supplied or canonical coordinates.

The patch does not modify:

- the custom best-conformer/current-reference reseeding between CREST iterations;
- COM or RMSD metadynamics mathematics;
- the separate `-TMD` scheduler budget;
- GFN-FF, EEQ, D3, fragment, optimizer, or external-rerank equations and semantics;
- input frame order or relative host–guest placements.

## Validation summary

The same source change was propagated to three private CREST 3.0.2 checkpoints:

- exact frozen-host caches;
- exact EEQ without the final MTD synchronization cleanup;
- exact EEQ with the final MTD synchronization cleanup.

Private validation established that all three builds:

- accepted five positional inputs with preoptimization enabled;
- retained the optimized first seed and preoptimized the four additional inputs;
- wrote numerically/byte-identical start ensembles at writer precision;
- constructed the expected `5 × 6 = 30` first-batch MTD jobs;
- preserved the single-input early-return behavior.

Negative tests covered changed atom count, changed element/order, inconsistent frozen-host coordinates, distorted guest topology, non-finite coordinates, and a reflected canonical target that cannot satisfy the proper-rotation postcondition.

No private structures, trajectories, chemical identities, or project directories are included in this publication.
