# zurand.so exports R_init_zurand and nothing else ($(C_VISIBILITY) in
# src/Makevars.in). Other packages reach the C API through the table
# registered with R_RegisterCCallable(), which needs no exported symbol;
# test-c-api.R covers that path.

test_that("the shared object exports exactly R_init_zurand", {
  skip_on_cran()
  skip_on_os("windows")
  # covr links the gcov runtime into the shared object, which exports its own
  # symbols.
  skip_if(nzchar(Sys.getenv("R_COVR")), "coverage build")
  skip_if(!nzchar(Sys.which("nm")), "nm not available")
  path <- getLoadedDLLs()[["zurand"]][["path"]]
  syms <- exported_symbols(path)
  skip_if(is.null(syms))
  # A coverage build links an instrumentation runtime (libgcov, LLVM's
  # profile runtime) that exports symbols of its own, not all of them
  # recognisable by name; that is not the shared object CRAN builds, and
  # every R-CMD-check leg runs this test on one that is.
  skip_if(any(grepl("gcov|llvm_prf|llvm_profile", syms)), "instrumented build")
  expect_identical(syms, "R_init_zurand")
})
