## Submission

This is a first submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

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
  Research), NumPy's ziggurat tables (NumPy Developers; the header is
  trimmed to the normal-double tables the package uses) and fdlibm's exp
  and log1p (Sun Microsystems).

## Test environments

GitHub Actions: macOS (arm64), Windows and Ubuntu (R release and
oldrel-1), Ubuntu arm64 (GCC), R-devel with clang and GCC 16 (r-hub
containers), plus ASan, UBSan, valgrind, rchk, LTO and gctorture; weekly
32-bit i386 and musl.
