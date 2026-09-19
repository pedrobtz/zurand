# Review of the Last Two Commits

Reviewed commits:

- `eff4160` - Speed up ziggurat normals with two-pass fill and derived wedge bounds
- `555591c` - Add `rng_threads()`, speed up uniform, fix wedge comment

## Findings

### 1. Low: the ziggurat shortcut is not fully conservative

`src/zurand.c:331` accepts or rejects without a margin on one side of the
fixed-point comparison, although `dev/generate-zig-bounds.R:98` acknowledges
that rounding of the `ki` values can move the curve across the chord.

For example, at layer 7 with `rabs == ki[7]`, the shortcut always accepts
because `R == L` and `YL < L`. The original `exp()` test rejects the top 188
possible 53-bit ordinates for that same point. The probability of observing
such a disagreement is extraordinarily small, but it contradicts the stated
invariant that the shortcut can only defer to the fallback and never disagree
with it.

Recommendation: add a small safety guard to both nominally margin-free
branches and send the boundary band through the `exp()` fallback. This should
have negligible performance cost.

### 2. Low: `rng_threads()` rejects values allowed by its documentation

`R/zurand.R:151` documents `threads` as any whole number of at least 1 and says
values above the machine thread count are allowed. `src/zurand.c:623` imposes an
undocumented upper bound of 1,048,576, while its error message only mentions
the lower bound.

Recommendation: either document the upper bound and include it in the error
message, or accept values through `INT_MAX`.

### 3. Low: the documented save/restore pattern does not restore an uncapped state

`C_rng_threads()` returns the current effective maximum from
`zurand_threads()`, rather than the raw configured cap. Consequently:

1. An initially uncapped package with an OpenMP maximum of 8 returns 8 from
   `old <- rng_threads(1)`.
2. `rng_threads(old)` stores an explicit cap of 8, rather than restoring the
   uncapped state.
3. If another package later raises the OpenMP maximum to 16, zurand remains
   capped at 8.

This makes the documented restoration pattern incomplete and weakens the
intended isolation from process-wide OpenMP configuration changes.

Recommendation: expose an explicit reset/no-cap operation or return a state
token that distinguishes an uncapped configuration from a cap equal to the
current OpenMP maximum.

## Correctness Assessment

The two-pass normal fill, uniform conversion, deferred uniform scaling, and
OpenMP fill plumbing otherwise appear correct. No data races or output-order
dependencies were found in the parallel paths.

Focused testing confirmed bit-identical uniform, normal, and integer results
across thread settings 1, 2, 3, 8, and an above-machine setting, for both
single-key and multi-key draws.

## Verification

- `dev/generate-zig-bounds.R` regenerated `src/zigbounds.h` byte-for-byte.
- `R CMD check --no-manual` completed successfully and all tests passed.
- The check reported two unrelated existing warnings:
  - placeholder license metadata in `DESCRIPTION`;
  - unused macOS OpenMP linker flags during compilation.
- The CRAN incoming check was not completed because its network-facing stage
  stalled.

## Verdict

The implementation is operationally sound and passes its tests, but the three
low-severity contract and boundary issues above should be addressed before the
new behavior is treated as formally exact.
