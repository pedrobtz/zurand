## Submission

This is a first submission.

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  New submission

* checking pragmas in C/C++ headers and code ... NOTE
  File which contains pragma(s) suppressing diagnostics: 'src/zurand.c'

  The pragma is `#pragma GCC diagnostic ignored "-Wunused-const-variable"`,
  pushed and popped around one `#include` of NumPy's ziggurat constant
  tables (src/numpyzig/ziggurat_constants.h, BSD 3-clause, credited in
  inst/COPYRIGHTS). The header is vendored verbatim and defines tables for
  several distributions, of which the package uses the normal ones; the
  pragma keeps the unused tables from producing warnings without editing
  the vendored file. No other diagnostic is suppressed.

## Notes for the reviewer

* Possibly misspelled words in DESCRIPTION, if flagged, are names:
  xoshiro, Philox and Threefry are the random number generators the
  package implements, and Blackman, Vigna, Salmon, Moraes, Dror, Marsaglia
  and Tsang are the authors of the cited papers.

* Multiple threads: the samplers use OpenMP above 32,768 values when R was
  built with it. The tests cap the package at two threads unless NOT_CRAN
  is set (tests/testthat/setup.R), and examples and the vignette draw only
  small vectors.
* Floating point: the package keeps compilers from fusing multiply-adds
  with an empty inline asm barrier rather than a compiler flag, so
  src/Makevars contains no -f flags. It computes the ziggurat's exp() and
  log1p() with its own port of fdlibm (src/zurand_fdlibm.h, Sun's notice
  preserved, credited in DESCRIPTION and inst/COPYRIGHTS) so output does not
  depend on the platform's libm.
* The package bundles three third-party components, all credited as `cph`
  in Authors@R and described in inst/COPYRIGHTS: Random123 (D. E. Shaw
  Research), NumPy's ziggurat tables (NumPy Developers) and fdlibm's exp
  and log1p (Sun Microsystems).

## Test environments

GitHub Actions: macOS (arm64), Windows and Ubuntu (R release and
oldrel-1), Ubuntu arm64 (GCC), R-devel with clang and GCC 16 (r-hub
containers), plus ASan, UBSan, valgrind, rchk, LTO and gctorture; weekly
32-bit i386 and musl.
