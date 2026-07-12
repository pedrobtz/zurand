# Public interface of `Random123/philox.h`

Summary of the API exposed by the vendored header ([src/Random123/philox.h](../src/Random123/philox.h), upstream: DEShawResearch/random123), written for readers who are comfortable with R but not with C.

## Background: what kind of RNG this is

R's built-in generator (Mersenne-Twister behind `runif()` etc.) is **stateful**: there is a hidden state (`.Random.seed`), every draw mutates it, and the *n*-th value can only be reached by generating the n−1 values before it, in order.

Philox is a **counter-based RNG (CBRNG)**. There is no state at all — it is a pure function:

```
out = philox(ctr, key)
```

- `ctr` (the *counter*) — "which position am I asking for?" Think of it as an index into an enormous pre-computed table of random numbers.
- `key` — "which table am I reading from?" Different keys give completely unrelated sequences.
- `out` — a block of high-quality random bits, fully determined by `(ctr, key)`.

Calling it twice with the same inputs gives the same output; changing even one bit of either input produces an unrelated-looking output. That single property is what makes random access, parallelism, and reproducibility trivial: value #1,000,000 costs the same as value #1, and two machines that agree on `(ctr, key)` agree on the value.

Philox is also a **bijection** (a perfect shuffle): as `ctr` runs over all possible counter values, `out` also hits every possible output exactly once — no collisions, no lost states. It is built from cryptographic ideas ("Product HI LO Xor": multiply, xor, repeat for several *rounds*) but is deliberately **not** cryptographically secure — it trades security margin for speed. It passes the full TestU01 BigCrush statistical battery.

## Words, widths, and variant names

C integers come in fixed sizes. Two matter here:

- `uint32_t` — an unsigned 32-bit integer, values 0 … 2³²−1 (≈ 4.3 billion)
- `uint64_t` — an unsigned 64-bit integer, values 0 … 2⁶⁴−1 (≈ 1.8 × 10¹⁹)

("unsigned" = no negative values; arithmetic wraps around modulo 2^width.)

A variant name like `philox4x64` reads as: the counter is **4 words of 64 bits** (a 256-bit counter), and one call returns 4 × 64 = **256 random bits**. The key is always half as many words as the counter.

| Variant | Counter | Key | Output per call | Min. safe rounds¹ | Availability |
|---|---|---|---|---|---|
| `philox2x32` | 2 × `uint32_t` | 1 × `uint32_t` | 64 bits | 6 | always |
| `philox4x32` | 4 × `uint32_t` | 2 × `uint32_t` | 128 bits | 8² | always |
| `philox2x64` | 2 × `uint64_t` | 1 × `uint64_t` | 128 bits | 6 | needs 64-bit multiply³ |
| `philox4x64` | 4 × `uint64_t` | 2 × `uint64_t` | 256 bits | 7 | needs 64-bit multiply³ |

¹ Fewest rounds with no known statistical flaws, per the upstream docs (2011). The default is 10 everywhere — a deliberate safety margin. A *round* is one pass of the multiply/xor scrambling step; more rounds = more thoroughly mixed, slightly slower.
² Philox4x32 with 7 rounds once showed suspicious p-values in very long statistical tests; 8+ is considered clean.
³ See [Portability of the 64-bit variants](#portability-of-the-64-bit-variants).

rngat uses **`philox4x64`** with the default 10 rounds ("Philox4x64-10").

## C API

`philox.h` is a *header-only* library: there is nothing to link against — `#include`-ing the file gives you everything, implemented as small functions the compiler pastes directly into your code (`static inline`, so they cost nothing to call).

For each variant `philoxNxW` the header defines the following. (Read `N` = 2 or 4, `W` = 32 or 64, e.g. `philox4x64_ctr_t`.)

### Types

C has no classes; the closest thing is a **struct** — a named bundle of fields. Counters and keys are structs holding a small fixed-size array named `v`:

```c
typedef struct r123arrayNxW  philoxNxW_ctr_t;   /* { uintW_t v[N];  } */
typedef struct r123arrayNhxW philoxNxW_key_t;   /* Nh = N/2 words     */
typedef struct r123arrayNhxW philoxNxW_ukey_t;  /* same thing, see below */
```

(`typedef` just introduces a shorter alias for a type; these array structs actually live in `array.h`, which `philox.h` includes.)

So a `philox4x64_ctr_t` is nothing more than "four `uint64_t`s in a row". You reach the words with `.v[i]` (0-indexed, like everything in C):

```c
philox4x64_ctr_t c = {{index, domain, purpose, 0}};  /* set all 4 words */
c.v[0]                                               /* read the first  */
```

The double braces are C initialization syntax: outer braces for the struct, inner for the array inside it.

`_key_t` vs `_ukey_t` ("user key"): for Philox they are the *same type*, and the conversion function below is the identity. The distinction only exists so all Random123 generators share one calling convention (for Threefry/AES the conversion does real work).

### Functions

```c
/* Run exactly R rounds (R <= 16, checked with an assertion). */
philoxNxW_ctr_t philoxNxW_R(unsigned int R, philoxNxW_ctr_t ctr, philoxNxW_key_t key);

/* Turn a "user key" into a key. For Philox: returns its input unchanged. */
philoxNxW_key_t philoxNxWkeyinit(philoxNxW_ukey_t uk);
```

Note that structs are passed and returned **by value** — the function gets a copy and hands back a fresh struct; nothing is modified in place. That is what makes the API stateless.

### Convenience macro and rounds constant

```c
philoxNxW(c, k)    /* == philoxNxW_R(philoxNxW_rounds, c, k) */
philoxNxW_rounds   /* named constant == 10 (the default rounds) */
```

A C **macro** is a compile-time text substitution — `philox4x64(c, k)` literally expands to `philox4x64_R(10, c, k)` before compilation. This is the form normal callers use.

## Compile-time configuration

C libraries are often configured with `#define NAME value` lines placed *before* the `#include` — the header checks "has the user already defined this?" and only supplies its default otherwise. Everything here is overridable that way:

- `PHILOXNxW_DEFAULT_ROUNDS` — rounds used by the convenience macro (default 10).
- Round-function multipliers: `PHILOX_M2x32_0`, `PHILOX_M4x32_0/_1`, `PHILOX_M2x64_0`, `PHILOX_M4x64_0/_1`.
- Key-schedule increments (the "Weyl constants", derived from the golden ratio and √3−1): `PHILOX_W32_0/_1`, `PHILOX_W64_0/_1`.

There is no reason to touch the multipliers/Weyl constants outside of research; changing them changes every output the generator ever produces.

## Portability of the 64-bit variants

The heart of each Philox round is a "mulhilo": multiply two W-bit numbers and keep **both** halves of the double-width product (a 64×64 multiply mathematically yields 128 bits). Standard C has no portable 128-bit integer, so the header picks an implementation per platform, in order of preference:

1. hand-written assembly / compiler intrinsics (MSVC, CUDA, OpenCL),
2. the GCC/Clang extension type `__uint128_t` (available on all 64-bit x86 and ARM — every platform CRAN currently builds for), or
3. a portable pure-C fallback that emulates the wide multiply with four 32-bit multiplies (`R123_USE_MULHILO64_C99`) — slower, but works anywhere.

If none of these is enabled, the `x64` variants simply don't exist (the flag `R123_USE_PHILOX_64BIT` ends up 0 and their types/functions are not compiled). **rngat force-enables the fallback** in [src/Makevars](../src/Makevars) (`-DR123_USE_MULHILO64_C99=1`) so `philox4x64` builds even on unusual 32-bit targets; on normal machines the faster `__uint128_t` path still wins because it is earlier in the preference order.

## C++ API (not used by rngat)

The header also contains a C++ layer, compiled only when included from C++. It wraps the same functions in small function-object classes:

```cpp
r123::PhiloxNxW_R<ROUNDS>   // rng(ctr, key) via operator()
r123::PhiloxNxW             // = PhiloxNxW_R<10>
```

These exist to plug into Random123's C++ adapter headers (`MicroURNG.hpp`, `uniform.hpp`, …), none of which are vendored into rngat.

## Minimal C example, line by line

```c
#include "Random123/philox.h"

philox4x64_key_t k = {{seed, 0}};                /* 128-bit key: word 0 = seed, word 1 = 0   */
philox4x64_ctr_t c = {{index, domain, purpose, 0}}; /* 256-bit counter, one meaning per word  */
philox4x64_ctr_t r = philox4x64(c, k);           /* 10 rounds; r is a fresh struct            */
/* r.v[0], r.v[1], r.v[2], r.v[3] are 4 independent uniform uint64s */
```

This is essentially the construction rngat uses (see `rngat_block()` in [src/rngat.c](../src/rngat.c)); the one refinement is that rngat packs **four** logical draws into each block instead of using only `r.v[0]` (see below). The R-facing consequences:

- `runif_at(key, index, domain)` evaluates the block at `{index >> 2, domain, 0, 0}`, picks word `index & 3`, and turns its top 53 bits into a double in (0, 1) — 53 because that is the precision of an R double's mantissa.
- `fold_in(key, id)` evaluates the block at `{id, 0, 1, 0}` — the third word (*purpose* = 1 instead of 0) guarantees derived keys can never overlap the bits used for draws — and uses `r.v[0]`, `r.v[1]` as the new 128-bit key.
- An `rngat_key` in R is just the two key words written out as 16 raw bytes (in a fixed byte order, so keys travel across machines).

## How rngat lays out the counter

Philox itself attaches no meaning to the counter — it is just 256 bits of input, and flipping *any* bit of it yields a statistically unrelated output block. rngat assigns one role to each 64-bit word:

| Word | Name | Set by | Meaning |
|---|---|---|---|
| `v[0]` | block index | user, as `index >> 2` (or `identity` for `fold_in`) | which 4-value block |
| `v[1]` | `domain` | user (`domain` argument, default 0) | which parallel stream at the same indices |
| `v[2]` | `purpose` | internal, never user-visible | what *kind* of question is being asked |
| `v[3]` | — | fixed 0 | reserved for future extensions |

`index` and `domain` are the user-facing coordinates: think of `(index, domain)` as addressing a cell in a 2⁶⁴ × 2⁶⁴ grid of pre-computed random values, any cell readable in any order.

### Four values per block

A `philox4x64` call produces 256 bits — four `uint64` words — but a single uniform or normal draw needs only 53 bits. Using just `r.v[0]` and discarding the other three words would waste ¾ of the generator's work. Instead rngat maps four consecutive logical positions onto one block:

```
value(index, domain) = word (index & 3) of philox( {index >> 2, domain, purpose, 0}, key )
```

So `index` is split into a **block number** (`index >> 2`, i.e. `index / 4`) placed in the counter, and a **word selector** (`index & 3`, i.e. `index % 4`) choosing which of the block's four words to return. Value #*i* is still an O(1) pure function of *i* — reproducible in any order, on any platform — but a contiguous run such as `runif_at(key, 1:1e6)` now costs about one Philox evaluation per **four** values. The C loop caches the most recently computed block, so neighbouring indices that share a block skip recomputation entirely; scattered indices simply fall back to one block each (never worse than before). This is the standard Random123 idiom (it is exactly what a sequential fill loop does — see randompack's `fill_philox`) and, because Philox's four output words are mutually independent, it changes none of the statistical guarantees.

### The `purpose` word: domain separation

`purpose` (constants `RNGAT_PURPOSE_*` in [src/rngat.c](../src/rngat.c)) partitions the counter space by *use*:

| Value | Constant | Used by |
|---|---|---|
| 0 | `RNGAT_PURPOSE_DRAW` | `bits_at()`, `runif_at()`, `rnorm_at()` |
| 1 | `RNGAT_PURPOSE_FOLD` | `fold_in()` |

**The bug it prevents.** Suppose `fold_in` simply evaluated Philox at `{id, 0, 0, 0}` — a counter in the same `purpose = 0` region that value draws use. Block `id` is exactly the block behind draws at indices `4*id … 4*id + 3` (domain 0), so `fold_in(key, id)` and those draws would read the *identical* output block:

```
fold_in(key, 2)        -> new key = {r.v[0], r.v[1]}   (first 128 bits of block 2)
bits_at(key, 8)        ->  value  =  r.v[0]            (word 0 of block 2 -- the same 64 bits!)
bits_at(key, 9)        ->  value  =  r.v[1]            (word 1 of block 2 -- also shared!)
```

The words of your derived key would literally *be* values you might also draw and print. Anything downstream of the derived key would be correlated with those draws — a subtle, hard-to-detect statistical defect. With `purpose = 1`, `fold_in` reads from a provably disjoint region of counter space: no output block is ever shared between key derivation and value draws, whatever `index`, `domain` or `identity` the user picks. (This mirrors what cryptographers call *domain separation*, the same reason protocols prefix hash inputs with distinct tags.)

**Extensibility.** The word also leaves room to grow: a future feature that needs its own dedicated bit-stream — say, a distribution that consumes more than 64 bits per variate, or an internal shuffling primitive — can claim `purpose = 2, 3, …` and is automatically independent of everything already defined. Because existing draws all live in the `purpose = 0` slice (and folds in `purpose = 1`), new purposes can be added without changing a single value the package produces today. The fixed `v[3] = 0` word is the same idea held in reserve: one entire spare 64-bit dimension, at zero cost.

**Rule for contributors:** any new use of `rngat_block()` that is not an ordinary value draw must use a fresh `purpose` value, never reuse 0 or 1.
