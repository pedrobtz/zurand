# Statistical audit

What `test-kat.R` and `test-golden.R` do not cover. Those prove the
generators really are Philox and Threefry, and that every CI platform
produces identical values. Neither says anything about whether the numbers
are *distributed* correctly, and the ziggurat -- the fixed-point wedge
shortcut and the tail in particular -- is this package's own code.

Reproduce with the `statistical-audit` workflow (on demand), or locally:

```sh
Rscript tools/statistical-audit.R --sampler=normal --engine=philox4x64 --gb=4 \
  | RNG_test stdin32
```

## Results

PractRand pre0.95, 256 MB per combination, x86_64, 2026-09-19.

| sampler | engine | result |
|---|---|---|
| normal | philox4x64 | no anomalies in 168 test results |
| normal | threefry4x64 | no anomalies in 168 test results |
| normal | xoshiro256pp | no anomalies in 180 test results (512 MB) |
| uniform | xoshiro256pp | no anomalies in 180 test results (512 MB) |
| normal | xoshiro256pp, after the fold fix | no anomalies in 168 test results (256 MB) |
| uniform | xoshiro256pp, after the fold fix | no anomalies in 168 test results (256 MB) |

The normal rows are the ones that matter: normals are mapped back through
`pnorm()` before testing, so a wrong ziggurat -- a mis-set wedge bracket, a
tail that returns the wrong branch -- shows up as non-uniformity. Both
engines are clean to 256 MB.

256 MB is a smoke test, not a verdict. PractRand's value is that it keeps
going: a generator that is clean at 256 MB and fails at 256 GB is
ordinary. The workflow takes a `gb` input for that reason, and a real run
before any release should be several hundred GB on the normal sampler for
both engines.

The xoshiro256pp rows matter for a second reason. That engine reseeds a
256-bit xoshiro state from Philox every 512 words, which is an unusual way
to drive xoshiro -- far more reseeding than its designers had in mind. If
that interacted badly with the recurrence, short-range correlation is
exactly what PractRand's lowest lengths would catch first. It does not, to
512 MB.

The last two rows exist because `tools/statistical-audit.R` derives one key
per block with `rng_fold()`, and a review found xoshiro keys were being
folded by the threefry path. Fixing that changed every block key the audit
uses for this engine, so the earlier rows no longer describe the stream the
tool emits; these do.

## Still to run

- uniform and bits, both engines, at the same depth
- a long run (>= 256 GB) on normal, which needs hours on a quiet machine
- dieharder as a second opinion, since batteries disagree about what is
  suspicious

## Reading the output

A stray WEAK among a hundred-plus tests per length is expected and means
nothing on its own. What matters is a FAIL, or a WEAK that returns at the
same test as the length doubles. The workflow summarises counts rather
than passing or failing the run, for that reason.
