# zurand (development version)

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
* `rng_threads()` caps OpenMP use; `rng_simd()` reports or disables the
  AVX2 path. Neither changes any value.

## C API

* `inst/include/zurand.h` lets packages that `LinkingTo: zurand` fill
  their own buffers, from their own threads, with exactly the values the R
  functions return.
