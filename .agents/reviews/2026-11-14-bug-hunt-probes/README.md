# Probes, bug hunt 2026-11-14

The trigger scripts that reproduced a finding, copied from the hunt's
scratch directory. `HUNT` stands for that directory, which held:

- `HUNT/lib`: zurand at `8a4c4eb` installed for macOS (serial, AVX2);
- `HUNT/linuxlib`: the same commit installed in
  `ghcr.io/r-hub/containers/ubuntu-release` (gcc 15, OpenMP, 4 threads),
  run as `podman run --rm -v HUNT:/w <image> Rscript /w/...`;
- `HUNT/srclinux`: `git archive 8a4c4eb`.

`coordinator/` holds the re-runs the report's statuses rest on.
`hunter-3/zc3` is the C API client used for #132 and #133.
