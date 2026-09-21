# Explicit charged GFN-FF fragments (CREST 3.0.2 fork)

This branch is based directly on `crest-3.0.2-final` and adds explicit formal
charges for the already-supported GFN-FF fragment selections.

Example:

```toml
[[calculation.level]]
method = "gfnff"
chrg = 0
uhf = 0
fragments = ["1-808", "809-938"]
fragment_charges = [-1, 1]
```

Rules:

- `fragment_charges` is optional. If absent, the validated 3.0.2-final behavior is unchanged.
- When present, there must be exactly one integer charge per explicit fragment.
- The fragment charges must sum exactly to `chrg`.
- Charged fragment selections must cover every atom exactly once and may not overlap.
- Charges are passed to GFN-FF as topology EEQ fragment constraints (`topo%qfrag`).
- The GFN-FF topology restart stores `qfrag`; process-isolated MTD also serializes
  the requested charge vector and validates it against the runtime topology.
- The process-MTD capsule format is bumped from V5 to V6 because the serialized
  calculator state gained the fragment-charge vector.

The original `crest-3.0.2-final` branch is not modified by this feature branch.
