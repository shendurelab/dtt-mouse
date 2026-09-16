# test_phylip_writer.R — standalone byte-identity + correctness check for the
# parallel dtt_phylip_writer.cpp (issue #40). Needs no project data: it builds
# synthetic SYMMETRIC matrices, so it validates the writer in isolation and runs
# anywhere (unlike test_full_distance_pipeline.R, which needs a real matrix).
#
# What it proves:
#   * the parallel writer's bytes are IDENTICAL to an independent, serial R-side
#     formatting of the same lower-triangular PHYLIP (the authoritative check);
#   * it crosses the internal CHUNK boundary (n > 4096) under real threads;
#   * NA entries are imputed to na_value and counted correctly;
#   * values round-trip at both ndigits=6 and the production ndigits=4.
#
# Run from the repo root (tree_building/):
#   Rscript 1_build_nj_backbone/tests/test_phylip_writer.R
suppressMessages({ library(Rcpp); library(RcppParallel) })
RcppParallel::setThreadOptions(numThreads = 4)   # exercise real concurrency
sourceCpp("1_build_nj_backbone/dtt_phylip_writer.cpp")

pass <- 0; fail <- 0
ok <- function(name, cond, extra = "") {
  if (isTRUE(cond)) { pass <<- pass + 1; cat(sprintf("PASS  %-40s %s\n", name, extra)) }
  else              { fail <<- fail + 1; cat(sprintf("FAIL  %-40s %s\n", name, extra)) }
}

# Independent reference: format the lower-triangular PHYLIP exactly as the C
# writer should, purely in R, to assert byte-identity.
ref_bytes <- function(M, labels, na_value, ndigits) {
  n <- nrow(M); fmt <- sprintf(" %%.%df", ndigits)
  lines <- character(n + 1L)
  lines[1] <- sprintf("%d", n)
  for (i in seq_len(n)) {
    vals <- ""
    if (i > 1) {
      v <- M[i, 1:(i - 1)]; v[is.na(v)] <- na_value
      vals <- paste0(sprintf(fmt, v), collapse = "")
    }
    lines[i + 1L] <- paste0(labels[i], vals)
  }
  paste0(paste0(lines, collapse = "\n"), "\n")
}

make_sym <- function(n, seed, na_frac = 0) {
  set.seed(seed)
  M <- matrix(runif(n * n, 0, 12), n, n)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]; diag(M) <- 0
  if (na_frac > 0) {
    k <- floor(na_frac * n * (n - 1) / 2); lt <- which(lower.tri(M))
    sel <- sample(lt, k)
    M[sel] <- NA; M <- t(M); M[sel] <- NA; M <- t(M)   # symmetric NA pattern
  }
  M
}

check_case <- function(tag, n, seed, na_frac, ndigits) {
  M <- make_sym(n, seed, na_frac)
  labels <- sprintf("cell_%05d", seq_len(n))
  na_value <- 3.14159
  plain <- tempfile(fileext = ".phy"); on.exit(unlink(plain), add = TRUE)

  n_imp <- write_phylip_lower(M, labels, na_value = na_value, path = plain, ndigits = ndigits)
  got <- readChar(plain, file.info(plain)$size, useBytes = TRUE)
  ok(sprintf("[%s] byte-identical to R reference", tag),
     identical(got, ref_bytes(M, labels, na_value, ndigits)))
  exp_imp <- sum(is.na(M[lower.tri(M)]))
  ok(sprintf("[%s] imputed count = %d", tag, exp_imp), n_imp == exp_imp)

  rows <- strsplit(trimws(readLines(plain)[2:(n + 1)]), "[[:space:]]+")
  back <- matrix(0, n, n)
  for (r in 2:n) back[r, 1:(r - 1)] <- as.numeric(rows[[r]][-1])
  back[upper.tri(back)] <- t(back)[upper.tri(back)]
  Mi <- M; Mi[is.na(Mi)] <- na_value
  tol <- 10^(-ndigits)
  ok(sprintf("[%s] round-trips within %.0e", tag, tol),
     max(abs(Mi[lower.tri(Mi)] - back[lower.tri(back)])) <= tol,
     sprintf("(n=%d, ndigits=%d)", n, ndigits))
}

check_case("small/6dig",      50,   1, 0.00, 6)
check_case("multichunk/6dig", 5000, 2, 0.00, 6)   # > CHUNK=4096 -> chunk loop + threads
check_case("with-NA/6dig",    400,  3, 0.05, 6)
check_case("prod/4dig",       5000, 4, 0.03, 4)   # production config (run_full_distance.R)

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(status = 1)
