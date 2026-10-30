# Helpers for test-abi.R, kept out of it so that the tests run shuffled.

exported_symbols <- function(path) {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Linux") {
    out <- system2("nm", c("-D", "--defined-only", shQuote(path)), stdout = TRUE)
  } else if (sys == "Darwin") {
    out <- system2("nm", c("-gU", shQuote(path)), stdout = TRUE)
  } else {
    return(NULL)
  }
  syms <- sub("^.* ", "", out)
  syms <- sub("^_(R_init)", "\\1", syms)
  # Linker-defined section markers are not code.
  setdiff(syms, c("", "_edata", "_end", "__bss_start", "_init", "_fini",
                  "__end__", "__bss_start__", "_bss_end__", "__bss_end__"))
}
