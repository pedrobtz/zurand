# Comprehensive RNG benchmark: zurand against the popular R generators.
#
# Run against the INSTALLED package (R CMD INSTALL, -O2); a
# pkgload::load_all() debug build understates zurand -- see CLAUDE.md.
#
#   R CMD INSTALL . && Rscript tools/benchmark.R
#
# Most competitors are single-threaded at the R level (base R, randompack,
# dqrng, sitmo, RcppZiggurat); rTRNG also parallelizes, via its parallelGrain
# argument. Each metric is shown twice: first with the threading generators
# (zurand, rTRNG) held to one thread -- algorithm vs algorithm -- then with
# threading allowed. The purely single-threaded packages are repeated in both
# tables unchanged, as a fixed yardstick.

library(zurand)

# Results are only comparable across machines if you know which machine
# produced them. zurand is tuned on Apple Silicon; an x86_64 run differs in
# cache line size, in whether Philox 64x64->128 is one instruction (x86 mulq)
# or two (arm64 mul + umulh), and in whether OpenMP is present at all.
platform <- function() {
  arch <- R.version$arch
  cpu <- tryCatch(
    if (Sys.info()[["sysname"]] == "Darwin")
      system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
    else sub("^model name\\s*:\\s*", "",
             grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
    error = function(e) NA_character_, warning = function(e) NA_character_)
  if (is.na(cpu[1]) || !nzchar(cpu[1]))
    cpu <- if (arch == "aarch64") "Apple Silicon" else arch
  cat("--------------------------------------------------------------\n")
  cat("CPU      : ", cpu[1], "\n", sep = "")
  cat("arch     : ", arch, " | cores: ", parallel::detectCores(), "\n", sep = "")
  cat("R        : ", R.version.string, "\n", sep = "")
  cat("zurand   : ", as.character(utils::packageVersion("zurand")), "\n", sep = "")
  cat("OpenMP   : ", if (rng_threads() > 1)
        paste0("yes, up to ", rng_threads(), " threads")
      else "NO -- single-threaded build, threaded tables below are meaningless",
      "\n", sep = "")
  cat("--------------------------------------------------------------\n")
}
platform()

required <- c("randompack", "dqrng", "sitmo", "RcppZiggurat", "rTRNG", "bench")
missing <- required[!vapply(required, requireNamespace, logical(1),
                            quietly = TRUE)]
if (length(missing))
  stop("install for the full benchmark: ", paste(missing, collapse = ", "))

suppressMessages({
  library(dqrng)
  library(sitmo)
  library(RcppZiggurat)
  library(rTRNG)
})

n <- 1e7L
seed <- 42L
TRNG_GRAIN <- 65536L  # rTRNG chunk size; parallelGrain = 0 runs it serially

# seed every generator so runs are reproducible
key <- rng_key(seed, engine = "philox4x64")
key_tf <- rng_key(seed, engine = "threefry4x64")
key_xo <- rng_key(seed, engine = "xoshiro256pp")
# randompack's default engine is its fastest, an 8-lane SIMD xoshiro256++;
# its philox engine is kept as the counter-based comparison.
rp <- randompack::randompack_rng()
rp_philox <- randompack::randompack_rng("philox")
set.seed(seed)
dqset.seed(seed)
zsetseed(seed)
TRNGseed(seed)

# --- pretty-printer: sort by median, show throughput and relative speed ------
show <- function(title, res) {
  res <- res[order(as.numeric(res$median)), ]
  med <- as.numeric(res$median)
  tab <- data.frame(
    method    = as.character(res$expression),
    median    = format(res$median),
    `M vals/s`= round(n / med / 1e6, 1),
    slowdown  = sprintf("%.2fx", med / min(med)),
    check.names = FALSE
  )
  cat("\n== ", title, "  (n = ", format(n, big.mark = ","), ") ==\n", sep = "")
  print(tab, row.names = FALSE)
  invisible(tab)
}

# env = parent.frame() so expressions see the caller's locals (e.g. the
# rTRNG grain), not this wrapper's frame.
mark <- function(...)
  bench::mark(..., check = FALSE, min_iterations = 20, env = parent.frame())

restore_threads <- rng_threads()
on.exit(rng_threads(restore_threads))

cat("zurand effective max threads:", rng_threads(), "\n")
cat("zurand SIMD path:", rng_simd(), "\n")
# `zurand` below is the default philox4x64 engine; the other two rows are
# the opt-in engines, and xoshiro256pp is the fast one.

# ============================ UNIFORM =======================================
uniform_bench <- function(trng_grain) mark(
  zurand_xoshiro  = rng_uniform(key_xo, n),
  zurand_threefry = rng_uniform(key_tf, n),
  zurand      = rng_uniform(key, n),
  randompack = rp$unif(len = n),
  randompack_philox = rp_philox$unif(len = n),
  dqrng      = dqrunif(n),
  base_R     = runif(n),
  sitmo      = runif_sitmo(n, 0, 1, seed),
  rTRNG      = runif_trng(n, parallelGrain = trng_grain)
)

rng_threads(1L)
show("UNIFORM, single thread", uniform_bench(0L))

rng_threads(0L)
show("UNIFORM, threading allowed (only zurand & rTRNG scale)",
     uniform_bench(TRNG_GRAIN))

# ============================ NORMAL ========================================
# base R rnorm() uses inversion (base R has no ziggurat for normals). The
# ziggurat spread comes from RcppZiggurat's classic variants:
#   zrnorm    - Leong et al. corrected ziggurat (its default)
#   zrnormMT  - Marsaglia & Tsang original
#   zrnormGSL - GNU Scientific Library
#   zrnormQL  - Gretl / QuantLib
# plus dqrng and randompack, which are ziggurat too.
normal_bench <- function(trng_grain) mark(
  zurand_xoshiro  = rng_normal(key_xo, n),
  zurand_threefry = rng_normal(key_tf, n),
  zurand          = rng_normal(key, n),
  randompack     = rp$normal(len = n),
  randompack_philox = rp_philox$normal(len = n),
  dqrng          = dqrnorm(n),
  base_inversion = rnorm(n),
  RcppZig_LZLLV  = zrnorm(n),
  RcppZig_MT     = zrnormMT(n),
  RcppZig_GSL    = zrnormGSL(n),
  RcppZig_QL     = zrnormQL(n),
  rTRNG          = rnorm_trng(n, parallelGrain = trng_grain)
)

rng_threads(1L)
show("NORMAL, single thread", normal_bench(0L))

rng_threads(0L)
show("NORMAL, threading allowed (only zurand & rTRNG scale)",
     normal_bench(TRNG_GRAIN))
