# CREST 3.0.2 host–guest workflow patch set

This directory publishes the exact validated changes developed for the CREST 3.0.2-v2 host–guest workflow.

## Baseline

- CREST repository commit: `cfdc301f759686b0fd66ced63b5ddbd6c693fa4f`
- Validated local source head: `092377fe68dd6cf493c98c377aca8861213d859c`
- GFN-FF is kept as the upstream CREST submodule; its three-file compatibility backport is applied as a patch.

## Apply

From a checkout of the baseline commit with submodules populated:

```bash
git checkout cfdc301f759686b0fd66ced63b5ddbd6c693fa4f
git submodule update --init --recursive
bash patches/crest-3.0.2-host-guest/apply.sh
```

The script applies:

1. `crest-3.0.2-host-guest.patch` to the top-level CREST source;
2. `gfnff-host-guest-fragments.patch` to the local `subprojects/gfnff` checkout.

## Implemented features

- separate `-TMD` dynamics/metadynamics thread budget;
- ordinary COM metadynamics bias and current-reference reseeding;
- multi-input NCI placement support;
- explicit host/guest fragmentation for the CREST 3.0.2 GFN-FF backend;
- staged external-rerank restart that resumes from the selected original CREST frame.

The GFN-FF backport changes only fragment assignment/topology construction. It does not alter GFN-FF parameters or force-field equations.

## Validation

The source was compiled and tested as local commit `092377f`. The staged 12-frame external-rerank restart regression passed and reproduced the expected output byte-for-byte.

Validated executable SHA-256:

```text
ae1d54d0aeb5e019690e94aa553817137121e1d082fc147d1745350f8454fb0e
```

Complete source/executable archive SHA-256:

```text
4c1c7b7e83020bc53c63686071b52fb78cbcdecc05ad8753c68c9e2f94b478f2
```
