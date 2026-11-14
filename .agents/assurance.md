# Assurance

How this package meets `.claude/shared/assurance.md`. CI checks this
file (r-actions `assurance.yml`), so keep the headings and table
columns as they are.

- **Statuses:** a cell is a fact, `gap #<issue>` or `n/a: <reason>`.

## Kinds

`compiled` `untrusted-input` `spec` `lifecycle`

## Required checks

`Where` names the workflow file under `.github/workflows/` that runs the
check. Add `(own)` when that workflow is the package's own rather than a
call to r-actions.

| Check         | Status | Where              |
| ------------- | ------ | ------------------ |
| sanitizers    | done   | native-checks.yaml |
| valgrind      | done   | native-checks.yaml |
| rchk          | done   | native-checks.yaml |
| gctorture     | done   | native-checks.yaml |
| analyzers     | done   | native-checks.yaml |
| alloc-failure | gap #108 |                    |
| mutation      | gap #109 |                    |
| fuzz          | gap #110 |                    |
| fuzz-reach    | gap #111 |                    |

## Boundaries

One row per boundary the package crosses. An oracle knows the answer:

- a conformance suite or test vectors;
- a reference implementation;
- a round trip;
- equivalent inputs compared;
- the bytes the other side received.

"No crash" is not an oracle.

<!-- markdownlint-disable MD013 -->
| Boundary | Oracle | Carried by |
| --- | --- | --- |
| Philox4x64-10 bit stream | Random123's published `kat_vectors`, and a from-spec Python Philox validated against all six of its philox4x64 rows | tests/testthat/test-kat.R (dev/kat/philox_ref.py) |
| Threefry4x64-13 bit stream | a from-spec Python Threefry | tests/testthat/test-kat.R (dev/kat/threefry_ref.py) |
| xoshiro256++ bit stream, Philox-seeded | a Python xoshiro256++ validated against Vigna's reference C | tests/testthat/test-kat.R (dev/kat/xoshiro_ref.py) |
| `rng_bits(bits = 32)` | equals the low half of the 64-bit stream | tests/testthat/test-kat.R |
| uniform transform (word to double, scaling, bounds more than `DBL_MAX` apart) | gap #112 | gap #112 |
| ziggurat normal: wedge and edge shortcut decisions | equal the `exp()` test computed in R from the defining tables | tests/testthat/test-shortcuts.R |
| ziggurat normal: values (layer, sign, wedge, tail) | gap #112 | gap #112 |
| McFarland normal: values (layers, edges, tail) | gap #112 | gap #112 |
| bounded integers (Lemire rejection) | gap #112 | gap #112 |
| fdlibm `exp()` and `log1p()` port | gap #113 | gap #113 |
| cross-platform bit stability | every CI leg reproduces golden values and whole-stream digests recorded on the reference platform | tests/testthat/test-golden.R, tests/testthat/test-digest.R |
| SIMD path against portable path | identical values (equivalent inputs) | tests/testthat/test-simd.R |
| threaded fill against serial fill | identical values under `rng_threads()` caps | tests/testthat/test-zurand.R |
| positional draws (`offset`) | equal the tail of a longer fill, on and off chunk boundaries | tests/testthat/test-zurand.R |
| C API (`LinkingTo: zurand`): fills, positional fills, streams | equal the R samplers, from worker threads, through a fixture package | tests/testthat/test-c-api.R |
| `rng_fold()`: typed data to key | integers and whole-number doubles fold identically; split points differ; snapshots | tests/testthat/test-zurand.R, tests/testthat/test-golden.R |
| `rng_fold()`: one string in two encodings | gap #114 | gap #114 |
| lazy vector save and restore | round trip equals the sampler; tampered states are refused | tests/testthat/test-lazy.R |
| key save and restore (stream attribute) | a saved key restores as its stream; malformed ones are refused | tests/testthat/test-key.R |
| statistical quality | PractRand through zurngtest's C callable | tests/testthat/test-zurngtest.R |
<!-- markdownlint-enable MD013 -->

## Objects

One row per object that outlives the call that made it. "Sequences" is
the test of its out-of-order use (abort then read, use after end, use
after free). The guard columns name the test, or say why the guard
cannot apply.

<!-- markdownlint-disable MD013 -->
| Object | Sequences | Fork | Unwind | Owner |
| --- | --- | --- | --- | --- |
| lazy vector (`rng_lazy_uniform()`, `rng_lazy_normal()`) | tests/testthat/test-lazy.R: read every way, write then save, read then save, materialise once, copies stay independent, interrupt mid-scan | n/a: holds no OS resource or thread, only R vectors | n/a: all memory is R's (RAWSXP state, `R_alloc`), so a longjmp leaks nothing; interrupt tested in tests/testthat/test-lazy.R | n/a: holds no borrowed pointer; the key is copied by value in `lazy_make()` |
| OpenMP thread pool (process-wide) | tests/testthat/test-zurand.R: `rng_threads()` caps and restores, values unchanged | gap #115 | n/a: parallel regions call no R API; interrupts are checked on the R thread | n/a: owned by the OpenMP runtime |
<!-- markdownlint-enable MD013 -->

## Limits

One row per limit the package enforces. "Straddled by" is the test of
`limit - 1`, `limit` and `limit + 1` on every path that reaches it.

<!-- markdownlint-disable MD013 -->
| Limit | Unit | Straddled by |
| --- | --- | --- |
| `offset + n` at most 2^53 | positions | tests/testthat/test-zurand.R ("offset + n is checked exactly at 2^53") |
| engine chunk (512 words; 4096 or 5120 for xoshiro256pp) | words | tests/testthat/test-golden.R, test-zurand.R, test-lazy.R, test-c-api.R |
| SIMD group (4 sub-chunks of 512) | values | tests/testthat/test-simd.R |
| C API positions up to 2^64 | positions | tests/testthat/test-c-api.R (starts at 2^64 - 2^16) |
| lazy index past 2^31 | elements | tests/testthat/test-lazy.R |
| key seed and `rng_fold()` doubles at most 2^53 - 1 | value | gap #116 |
| `rng_integer()` bounds in R's non-missing integer range | value | gap #116 |
| OpenMP threshold (`ZURAND_OMP_MIN_VALUES`, 32768) | values | gap #117 |
| huge-page advice (`ZURAND_THP_MIN_BYTES`, 4 MiB) | bytes | gap #117 |
| lazy access paths (`ZURAND_LAZY_RUN` 16, `ZURAND_LAZY_BIG_REGION` 2^15, `ZURAND_LAZY_COMPARE` 2^16) | elements | gap #118 |
<!-- markdownlint-enable MD013 -->
