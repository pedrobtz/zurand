# zurand: Fast Stateless Random Number Generation

Fast, stateless, deterministic random number generation. Random values
are pure functions of immutable keys and sampler arguments, so sampling
never reads or mutates '.Random.seed', and output is bit-identical
across platforms and thread counts. The default engine is xoshiro256++
(Blackman and Vigna, 2021)
[doi:10.1145/3460772](https://doi.org/10.1145/3460772) , reseeded from
the Philox counter-based generator of the 'Random123' library (Salmon,
Moraes, Dror and Shaw, 2011)
[doi:10.1145/2063384.2063405](https://doi.org/10.1145/2063384.2063405) ;
Philox and Threefry are also available directly. Normal deviates use the
ziggurat method (Marsaglia and Tsang, 2000)
[doi:10.18637/jss.v005.i08](https://doi.org/10.18637/jss.v005.i08) .
Independent substreams are requested as key vectors with rng_key(seed,
n) or derived explicitly with rng_fold().

## See also

Useful links:

- <https://pedrobtz.github.io/zurand/>

- <https://github.com/pedrobtz/zurand>

- Report bugs at <https://github.com/pedrobtz/zurand/issues>

## Author

**Maintainer**: Pedro Baltazar <pedrobtz@gmail.com> \[copyright holder\]

Authors:

- Pedro Baltazar <pedrobtz@gmail.com> \[copyright holder\]

Other contributors:

- D. E. Shaw Research (Random123 library, src/Random123/) \[copyright
  holder\]

- NumPy Developers (Ziggurat constant tables, src/numpyzig/) \[copyright
  holder\]

- Sun Microsystems, Inc. (fdlibm exp() and log1p(), src/zurand_fdlibm.h)
  \[copyright holder\]

- Christopher D. McFarland (modified ziggurat (fast_prng), from which
  tools/generate-mcfarland-tables.py and rng_normal(method =
  'mcfarland') derive) \[copyright holder\]
