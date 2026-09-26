# SIMD microbenchmarks

The measurements behind Phase 5 of `dev/roadmap.md`. Standalone C, no R.

```sh
cd src
clang -O2 -falign-functions=64 -DR123_USE_MULHILO64_C99=1 -I. \
      -o /tmp/lanes  ../dev/simd/lanes.c  && /tmp/lanes  10000000   # baseline
clang -O2 -falign-functions=64 -mavx2 -DR123_USE_MULHILO64_C99=1 -I. \
      -o /tmp/lanes  ../dev/simd/lanes.c  && /tmp/lanes  10000000   # AVX2
clang -O2 -falign-functions=64 -mavx2 -DR123_USE_MULHILO64_C99=1 -I. \
      -o /tmp/xchunk ../dev/simd/xchunk.c && /tmp/xchunk 10000000
```

`lanes.c` compares word sources inside zurand's two-pass fill: the shipped
one-lane recurrence against 2, 4 and 8 lanes interleaved within a chunk,
plus a hand-written AVX2 4-lane version, plus a control that pays eight
Philox seeds per chunk to isolate seeding cost from register pressure.

`xchunk.c` is the one that mattered: four consecutive chunks in four AVX2
lanes, each running the shipped sequential recurrence, with a 4x4
transpose so every chunk gets contiguous stores. It asserts word-for-word
equality with the scalar engine before timing anything.

Parse the output on the field before `M/s`, not on `$3` -- labels contain
spaces, and both earlier attempts to summarise these runs with `awk`
silently read the wrong column.

## Second pass, x86 (2026-09-26)

The harnesses behind the roadmap's "x86 analysis, second pass". All
self-contained C except the two R scripts; build lines are in each header.

- `fuse_normal.c`: where a normal fill's time goes and whether fusing the
  AVX2 generator with the ziggurat transform helps (it does not), plus
  non-temporal stores on fresh memory (slower). Reused buffer and fresh
  `malloc()` per repetition, as R allocates.
- `pagefault.c`: what a fresh 80 MB vector costs on this kernel, with and
  without huge pages (`MADV_HUGEPAGE`, `MADV_NOHUGEPAGE`,
  `MADV_POPULATE_WRITE`). The 4 KiB-page row is what a `madvise`-mode
  system pays without the hint that #42 added.
- `avx512.c`: 8-lane AVX-512 xoshiro with `vprolq` and an 8x8 transpose,
  asserted word-equal to the AVX2 path, round robin against it. Runs only
  where the CPU has avx512f and avx512dq.
- `zig_slow.c`: the slow path taken apart -- how often a rejected draw
  reaches `exp()` (about 6% of wedge tests), and `zurand_exp`/`log1p`
  against libm.
- `norm_ab.R`: normal-only same-machine A/B of two installed builds
  (`lib-a`, `lib-b`), one thread, ten alternating rounds.
- `ab_thp.R`: the huge-page experiment at the R level, one build run with
  `ZURAND_THP=off|no|huge` (the switch existed only on that branch).
