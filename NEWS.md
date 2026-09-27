# zurand 0.1.0

First release. Stateless random numbers: every value is a pure function of
an immutable key and the draw's position, so output is bit-identical
across platforms, thread counts, SIMD paths and call order.

## The stream (frozen as `stream-1`)

* Three engines. `xoshiro256pp`, the default, is xoshiro256++ reseeded
  every 512 values from Philox4x64-10; `philox4x64` and `threefry4x64` are
  the Random123 counter-based generators, available directly.
* Every key records its engine and a stream version (`attr(key, "stream")`,
  currently 1), so a stream can be corrected after release without
  changing what existing keys produce.
* Normal draws use a ziggurat on NumPy's tables with zurand's own port of
  fdlibm's `exp()` and `log1p()`, and no fused multiply-add on any
  compiler, so they are identical on Linux (glibc and musl), macOS and
  Windows, on x86_64 and aarch64. 32-bit x87 is outside that guarantee.
* `xoshiro256pp` passed a PractRand audit clean to 1 TB, including
  interleaved streams of sibling and folded keys.

## Functions

* `rng_key()`, `rng_key_from_r()` and `rng_fold()` create and derive keys.
* `rng_uniform()`, `rng_normal()`, `rng_integer()` and `rng_bits()` sample;
  with a vector of keys they return one column per key.
* `rng_uniform()` and `rng_normal()` take `offset`: `rng_normal(key, n,
  offset = k)` returns positions `k` to `k + n - 1`, exactly
  `rng_normal(key, k + n)[(k + 1):(k + n)]`, so a long simulation can draw
  its next batch, or resume from a checkpoint, without the earlier values.
* `rng_lazy_uniform()` and `rng_lazy_normal()` return lazy vectors (ALTREP):
  created instantly, with each value computed only when read and exactly
  equal to the matching sampler's. `sum()` over 1e9 lazy normals runs in
  the memory R itself uses; saving one stores only the key and arguments.
  Arithmetic and C code that asks for the data pointer fill the vector
  once. For scattered reads, key it with the `philox4x64` engine.
* `rng_normal(method = "mcfarland")` samples with McFarland's (2016)
  modified ziggurat, whose common case needs no table comparison. Large
  fills run 1.15-1.25x faster than the default method on x86_64 and
  1.05-1.10x on arm64. It is a separate, equally reproducible stream; the
  default method is unchanged.
* `rng_threads()` caps OpenMP use; `rng_simd()` reports or disables the
  AVX2 path. Neither changes any value.

## Performance

* The `xoshiro256pp` engine vectorises its generator on x86_64 (AVX2,
  chosen at run time) and on arm64 (NEON), and converts uniforms in
  registers without a second pass. Output is identical with the vector path
  on or off, which `rng_simd()` can switch to check.
* The ziggurat's rare slow path is compiled as plain out-of-line code
  (fdlibm `exp()`/`log1p()` as calls, no `cold` attribute): 7-9% faster
  normals on a Xeon 6973P and 15-19% with Apple clang, same values.
* Large fills use OpenMP threads where R was built with it; `rng_threads()`
  caps them, and the values do not depend on the thread count. On macOS,
  where R does not set OpenMP flags, `configure` links the OpenMP runtime
  that R ships, as CRAN's macOS binaries of data.table do, and falls back
  to one thread when the OpenMP headers are missing.
* On Linux kernels with transparent huge pages in `madvise` mode, output
  vectors of 4 MiB and more are filled on 2 MiB pages, which about halves
  the time of large fills there. `options(zurand.hugepages = FALSE)` turns
  the request off.

## C API

* `inst/include/zurand.h` lets packages that `LinkingTo: zurand` fill
  their own buffers, from their own threads, with exactly the values the R
  functions return.
* API version 2 adds `fill_normal_method()`, which takes the normal method
  as `ZURAND_NORMAL_ZIGGURAT` or `ZURAND_NORMAL_MCFARLAND`, the C side of
  `rng_normal(method =)`.
* API version 3 draws from any position, in two ways. `fill_uniform_at()`
  and `fill_normal_at()` fill positions `start` to `start + n - 1`, so the
  caller pulls the next piece whenever it wants it. `stream_uniform()` and
  `stream_normal()` run the loop instead: they pour positions `0` to
  `total - 1` through a small buffer the caller owns, calling a consumer
  function on each chunk in order, so a long simulation consumes its
  numbers from cache without the whole stream ever being in memory. Either
  way the values are those of `rng_uniform()` and `rng_normal()`.
