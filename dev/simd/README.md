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
