#!/bin/sh
# Syntax-check the OpenMP pragmas on a machine that has no OpenMP.
#
# Apple clang ships no libomp, so on this Mac every #ifdef ZURAND_OPENMP
# block is compiled out and its pragmas are never parsed. That is how a
# variable used inside a default(none) region without a data-sharing clause
# reached CI: it built clean locally and failed on every Linux leg.
#
# dev/ompstub/omp.h declares just enough for -fsyntax-only to parse the
# pragmas. Nothing is linked and no code is produced.
#
#   sh dev/check-openmp.sh
set -e
RINC=$(Rscript -e 'cat(R.home("include"))')
cd "$(dirname "$0")/../src"
clang -fsyntax-only -Xclang -fopenmp \
      -I../dev/ompstub -I"$RINC" -I. -I../inst/include \
      -DNDEBUG -DR123_USE_MULHILO64_C99=1 zurand.c
echo "OpenMP pragmas OK (default(none) clauses complete)"
