#!/usr/bin/env Rscript
# test_shim_correctness.R -- FAST, clean correctness gate for the in-memory NJ shim.
# Fixes the oracle's flaw: compare the shim to an independent canonical Saitou-Nei
# NJ on FULL-PRECISION (unrounded) distance matrices, which are effectively tie-free
# (Qmargin >> 0 => the NJ tree is genuinely unique). Both the reference and the shim
# consume the SAME in-memory matrix (no text round-trip), so if the shim's NJ is
# correct the bipartition sets must be identical. Also catches crashes (tryCatch)
# and verifies no tips are dropped. Small n => completes in seconds.
suppressWarnings(suppressMessages(library(Rcpp)))
S <- "1_build_nj_backbone"; D <- Sys.getenv("DTT_DIR", "decenttree")
# OpenMP: Linux gcc/clang take -fopenmp directly. Apple clang / Homebrew's
# libomp don't resolve decenttree's OpenMP usage cleanly, so Darwin skips
# OpenMP entirely and compiles single-threaded instead (see run_full_distance.R).
if (Sys.info()[["sysname"]] == "Darwin") {
  omp_cxxflags <- ""
  omp_libs      <- ""
} else {
  omp_cxxflags <- "-fopenmp"
  omp_libs      <- "-fopenmp"
}
Sys.setenv(PKG_CXXFLAGS = paste0("-I", normalizePath(D),
                                 " -I", normalizePath(file.path(D, "build")), " -O2 ", omp_cxxflags))
Sys.setenv(PKG_LIBS = omp_libs)
sourceCpp(file.path(S, "nj_inmemory.cpp"))
sourceCpp(file.path(S, "dtt_phylip_writer.cpp"))   # CLI-path input, for the consistency test
DECENTTREE <- "decenttree/build/decenttree"

# --- independent canonical NJ (Saitou-Nei), returns clades + smallest Q-margin ---
nj_canonical <- function(D, tips) {
  active <- as.list(tips); Dm <- D; clades <- list(); min_margin <- Inf
  while (length(active) > 3) {
    r <- length(active); Sr <- rowSums(Dm)
    Q <- (r - 2) * Dm - outer(Sr, Sr, `+`); diag(Q) <- Inf; Q[lower.tri(Q)] <- Inf
    qv <- Q[is.finite(Q)]; if (length(qv) >= 2) { sq <- sort(qv); min_margin <- min(min_margin, sq[2]-sq[1]) }
    ij <- which(Q == min(Q), arr.ind = TRUE)[1, ]; i <- ij[1]; j <- ij[2]
    u <- c(active[[i]], active[[j]]); clades[[length(clades)+1]] <- u
    duk <- 0.5 * (Dm[i, ] + Dm[j, ] - Dm[i, j]); keep <- setdiff(seq_len(r), c(i, j))
    newD <- Dm[keep, keep, drop = FALSE]
    Dm <- rbind(cbind(newD, duk[keep]), c(duk[keep], 0)); active <- c(active[keep], list(u))
  }
  list(clades = clades, min_margin = min_margin)
}
newick_clades <- function(nwk) {
  s <- gsub(":[0-9eE.+-]+", "", gsub("[[:space:];]", "", nwk))
  stack <- list(); clades <- list(); token <- ""
  for (ch in strsplit(s, "")[[1]]) {
    if (ch == "(") stack[[length(stack)+1]] <- character(0)
    else if (ch == "," || ch == ")") {
      if (nchar(token) > 0) { L <- length(stack); stack[[L]] <- c(stack[[L]], token); token <- "" }
      if (ch == ")") { grp <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
        clades[[length(clades)+1]] <- grp
        if (length(stack) > 0) { L <- length(stack); stack[[L]] <- c(stack[[L]], grp) } }
    } else token <- paste0(token, ch)
  }
  clades
}
biparts <- function(clades, all_tips) {
  ref <- all_tips[1]; n <- length(all_tips); keys <- character(0)
  for (c in clades) { side <- if (ref %in% c) setdiff(all_tips, c) else c
    if (length(side) >= 2 && length(side) <= n - 2) keys <- c(keys, paste(sort(side), collapse = "|")) }
  unique(keys)
}

cat("=== shim NJ vs independent canonical NJ, FULL-PRECISION (tie-free) matrices ===\n")
allpass <- TRUE; ncase <- 0
for (cs in list(c(10,1),c(15,2),c(20,3),c(25,4),c(30,5),c(40,6),c(50,7),c(60,8),c(80,9),c(100,10))) {
  n <- cs[1]; seed <- cs[2]; set.seed(seed)
  M <- matrix(0, n, n); ut <- upper.tri(M)
  M[ut] <- runif(sum(ut), 0.01, 1.0); M <- M + t(M)   # unrounded -> effectively no Q ties
  lab <- sprintf("t%03d", seq_len(n)); dimnames(M) <- list(lab, lab)
  ref <- nj_canonical(M, lab); refbp <- biparts(ref$clades, lab)
  out <- tempfile(fileext = ".nwk")
  res <- tryCatch({ decenttree_nj_inmem(M, lab, na_value = 0.5, nthreads = 1,
                                        precision = 6, out_path = out, ndigits = 0); "ok" },
                  error = function(e) paste("CRASH:", conditionMessage(e)))
  if (!identical(res, "ok")) { cat(sprintf("  n=%3d seed=%d  %s -> FAIL\n", n, seed, res)); allpass <- FALSE; next }
  nwk  <- trimws(readChar(out, file.info(out)$size, useBytes = TRUE))
  sbp  <- biparts(newick_clades(nwk), lab)
  rf    <- length(union(refbp, sbp)) - length(intersect(refbp, sbp))
  ntips <- length(gregexpr("[(,][^,():]+:", nwk, perl = TRUE)[[1]])
  tie_free <- ref$min_margin > 1e-9
  ok <- (rf == 0 && ntips == n && length(refbp) == n - 3 && length(sbp) == n - 3)
  allpass <- allpass && ok; ncase <- ncase + 1
  cat(sprintf("  n=%3d seed=%d  tips=%d/%d  bp=%d/%d  RF=%d  Qmargin=%.2e%s  -> %s\n",
              n, seed, ntips, n, length(refbp), length(sbp), rf, ref$min_margin,
              if (tie_free) "" else " (TIED!)", if (ok) "PASS" else "FAIL"))
  unlink(out)
}
cat(sprintf("\nCORRECTNESS: %s  (%d cases)\n", if (allpass) "PASS" else "FAIL", ncase))

# --- NA imputation: B2 has NA pairs. Shim must impute NaN->med and match both the
#     canonical NJ on the pre-imputed matrix AND the shim fed a pre-imputed matrix. ---
cat("\n=== NA (R NA_real_) imputation correctness ===\n")
na_pass <- TRUE
for (cs in list(c(30,201,0.1), c(60,202,0.15), c(90,203,0.2))) {
  n <- cs[1]; seed <- cs[2]; na_frac <- cs[3]; set.seed(seed)
  M <- matrix(0, n, n); ut <- upper.tri(M); M[ut] <- runif(sum(ut), 0.01, 1.0); M <- M + t(M)
  lab <- sprintf("t%03d", seq_len(n)); dimnames(M) <- list(lab, lab)
  med <- median(M[ut])
  Mna <- M; k <- floor(sum(ut) * na_frac); idx <- sample(which(ut), k)
  Mna[idx] <- NA_real_; Mt <- t(Mna); Mna[t(ut)] <- Mt[t(ut)]   # mirror NA to lower tri
  Mref <- Mna; Mref[is.na(Mref)] <- med                          # R-preimputed reference
  o1 <- tempfile(fileext=".nwk"); o2 <- tempfile(fileext=".nwk")
  decenttree_nj_inmem(Mna,  lab, med, 1, 6, o1, ndigits = 0)     # shim imputes NaN->med
  decenttree_nj_inmem(Mref, lab, 0.5, 1, 6, o2, ndigits = 0)    # nothing to impute
  ref  <- nj_canonical(Mref, lab)                                # canonical NJ on imputed M
  t1 <- trimws(readChar(o1, file.info(o1)$size, TRUE)); t2 <- trimws(readChar(o2, file.info(o2)$size, TRUE))
  rf <- length(union(biparts(newick_clades(t1),lab), biparts(ref$clades,lab))) -
        length(intersect(biparts(newick_clades(t1),lab), biparts(ref$clades,lab)))
  ok <- identical(t1, t2) && rf == 0
  na_pass <- na_pass && ok
  cat(sprintf("  n=%3d na=%2.0f%%  shim-impute==R-preimpute:%s  RF-vs-canonical=%d -> %s\n",
              n, 100*na_frac, identical(t1,t2), rf, if (ok) "PASS" else "FAIL"))
  unlink(c(o1, o2))
}

# --- self-consistency at -nt 1 (deterministic) ---
cat("\n=== self-consistency, -nt 1, byte-identical ===\n")
self_pass <- TRUE
for (cs in list(c(100,301), c(300,302))) {
  n <- cs[1]; seed <- cs[2]; set.seed(seed)
  M <- matrix(0, n, n); ut <- upper.tri(M); M[ut] <- runif(sum(ut), 0.01, 1.0); M <- M + t(M)
  lab <- sprintf("t%03d", seq_len(n)); dimnames(M) <- list(lab, lab)
  o1 <- tempfile(fileext=".nwk"); o2 <- tempfile(fileext=".nwk")
  decenttree_nj_inmem(M, lab, 0.5, 1, 6, o1, 0); decenttree_nj_inmem(M, lab, 0.5, 1, 6, o2, 0)
  ok <- identical(readLines(o1), readLines(o2)); self_pass <- self_pass && ok
  cat(sprintf("  n=%3d  run1==run2: %s -> %s\n", n, ok, if (ok) "PASS" else "FAIL"))
  unlink(c(o1, o2))
}

# --- production consistency (adversary finding #1): on real-shaped p/q distances
#     (q<=11), the shim at ndigits=4 must reproduce the CLI's ndigits=4 tree. Full
#     precision (ndigits=0) would flip ties and diverge (RF up to ~44) -- that is the
#     bug this check guards against. ---
cat("\n=== production consistency: shim(ndigits=4) vs CLI(ndigits=4), p/q distances ===\n")
cons_pass <- TRUE
for (cs in list(c(40,401), c(80,402), c(150,403))) {
  n <- cs[1]; seed <- cs[2]; set.seed(seed)
  M <- matrix(0, n, n); ut <- upper.tri(M)
  q <- sample(2:11, sum(ut), replace = TRUE); p <- vapply(q, function(qq) sample.int(qq, 1), integer(1))
  M[ut] <- p / q; M <- M + t(M); lab <- sprintf("t%03d", seq_len(n)); dimnames(M) <- list(lab, lab)
  med <- round(median(M[ut]), 4)
  phy <- tempfile(fileext = ".phy"); write_phylip_lower(M, lab, na_value = med, path = phy, ndigits = 4)
  system2("gzip", c("-f", shQuote(phy)))
  cli <- tempfile(fileext = ".nwk")
  system2(DECENTTREE, c("-in", shQuote(paste0(phy, ".gz")), "-t", "NJ-R", "-nt", 1, "-f", 6,
                        "-out", shQuote(cli), "-no-banner"), stdout = FALSE, stderr = FALSE)
  shim <- tempfile(fileext = ".nwk"); decenttree_nj_inmem(M, lab, med, 1, 6, shim, ndigits = 4)
  cbp <- biparts(newick_clades(trimws(readChar(cli,  file.info(cli)$size,  TRUE))), lab)
  sbp <- biparts(newick_clades(trimws(readChar(shim, file.info(shim)$size, TRUE))), lab)
  rf  <- length(union(cbp, sbp)) - length(intersect(cbp, sbp))
  ok  <- (rf == 0); cons_pass <- cons_pass && ok
  cat(sprintf("  n=%3d  RF(shim_4dp, CLI_4dp)=%d/%d  -> %s\n", n, rf, n - 3,
              if (ok) "PASS (reproduces CLI)" else sprintf("FAIL (%d splits differ)", rf)))
  unlink(c(paste0(phy, ".gz"), cli, shim))
}

cat(sprintf("\n=== VERDICT: correctness=%s  NA=%s  self-consistent=%s  CLI-consistency=%s ===\n",
            if (allpass) "PASS" else "FAIL", if (na_pass) "PASS" else "FAIL",
            if (self_pass) "PASS" else "FAIL", if (cons_pass) "PASS" else "FAIL"))
if (!(allpass && na_pass && self_pass && cons_pass)) quit(status = 1)
cat("ALL SHIM GATES PASS\n")
