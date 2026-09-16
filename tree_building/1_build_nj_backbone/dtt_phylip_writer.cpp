// dtt_phylip_writer.cpp
//
// Stream an in-memory symmetric distance matrix to a LOWER-TRIANGULAR PHYLIP
// distance file that DecentTree reads directly (row i carries its i entries
// for columns 1..i-1, no diagonal -- half the size of a square matrix, and
// DecentTree accepts this triangular form as-is, so we never materialise the
// square matrix on disk). NA entries are imputed inline with `na_value` (the
// global median distance) -- DecentTree itself would silently turn NA into 0
// ("identical"), which we must avoid.
//
// Writes PLAIN text via buffered stdio -- no compression here. Compression is
// a separate, later step done by piping through pigz (multi-threaded gzip),
// which parallelises across cores instead of the single-threaded deflate that
// zlib's gzwrite() used to do inline; see run_full_distance.R.
//
// Returns the number of NA entries that were imputed (lower triangle only), so
// the caller can sanity-check it against the expected NA rate.
//
// PERFORMANCE (issue #40 -- "slow distance matrix write"). At real scale the
// lower triangle is ~n^2/2 values (~6.7e10 at n=367k), and the naive writer
// was BOTH single-threaded (one snprintf per value on one core, ~hours) AND
// cache-hostile: it read m[i + j*n] with i fixed and j varying, striding the
// column-major matrix by n*8 bytes per value, so almost every read was a cache
// / TLB miss on a >1 TB matrix. Two fixes here, both provably lossless:
//   1. Read m[j + i*n] instead of m[i + j*n]. The matrix is SYMMETRIC and the
//      kernel writes out[i,j] and out[j,i] from the SAME computed double
//      (see dtt_distance_scale.cpp), so this is bit-identical -- but now the
//      j-loop for a fixed row i walks column i CONTIGUOUSLY in memory.
//   2. Format the (independent) rows in parallel with RcppParallel into
//      per-row buffers, then fwrite them in order. The bytes are identical to
//      the old serial writer (same snprintf, same format string, same order);
//      only the formatting is spread across cores. Rows are processed in
//      chunks so the transient buffer memory stays bounded (~CHUNK * longest
//      row), negligible next to the matrix itself.

#include <Rcpp.h>
#include <RcppParallel.h>
#include <cstdio>
#include <charconv>
#include <string>
#include <vector>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
using namespace Rcpp;

// Floating-point std::to_chars needs libstdc++ >= GCC 11 or libc++ >= 14; the
// feature-test macro is only defined once that support is complete. If the
// build toolchain is older we fall back to snprintf -- correct everywhere, and
// on Linux/glibc it still parallelises fine (the per-call locale lock that
// serialises it is macOS/BSD-libc specific, not a concern on the run node).
#if defined(__cpp_lib_to_chars) && (__cpp_lib_to_chars >= 201611L)
#  define DTT_HAVE_TO_CHARS 1
#else
#  define DTT_HAVE_TO_CHARS 0
#endif

// Format v as " <fixed-ndigits>" into buf (>= 64 bytes); return bytes written.
// fmt must be the matching " %.<ndigits>f" string (used only in the fallback).
static inline int dtt_format_value(char* buf, double v, int ndigits,
                                   const char* fmt) {
#if DTT_HAVE_TO_CHARS
  buf[0] = ' ';                           // leading-space value separator
  auto r = std::to_chars(buf + 1, buf + 64, v, std::chars_format::fixed, ndigits);
  (void)fmt;
  return (int)(r.ptr - buf);
#else
  (void)ndigits;
  return std::snprintf(buf, 64, fmt, v);
#endif
}

// Formats a contiguous block of lower-triangular rows into per-row string
// buffers. Reads only the raw double* (thread-safe) and plain std::strings --
// never touches R objects from a worker thread. Uses std::to_chars (C++17):
// it is locale-independent and lock-free, unlike snprintf("%f") whose per-call
// locale lock serialised the old parallel attempt (and produces bit-identical
// fixed-notation output to "%.<ndigits>f" for these non-negative distances).
struct RowFormatter : public RcppParallel::Worker {
  const double* m;                       // n*n column-major distance matrix
  const std::vector<std::string>* labs;  // taxon labels (extracted in R thread)
  const int n;
  const double na_value;
  const int ndigits;
  const char* fmt;                       // " %.<ndigits>f" (snprintf fallback)
  const int base;                        // absolute row index of bufs[0]
  std::vector<std::string>* bufs;        // output: bufs[i-base] = row i's line
  std::vector<long long>* imp;           // output: imp[i-base] = NAs in row i

  RowFormatter(const double* m_, const std::vector<std::string>* labs_, int n_,
               double na_value_, int ndigits_, const char* fmt_, int base_,
               std::vector<std::string>* bufs_, std::vector<long long>* imp_)
    : m(m_), labs(labs_), n(n_), na_value(na_value_), ndigits(ndigits_),
      fmt(fmt_), base(base_), bufs(bufs_), imp(imp_) {}

  void operator()(std::size_t begin, std::size_t end) {
    char buf[64];
    for (std::size_t i = begin; i < end; ++i) {
      std::string& row = (*bufs)[i - base];
      row.clear();
      row.reserve((*labs)[i].size() + i * (std::size_t)(ndigits + 8) + 1);
      row.append((*labs)[i]);
      long long local_imp = 0;
      // Column i, rows 0..i-1 -- contiguous in memory (symmetry: M(i,j)==M(j,i)).
      const double* col = m + (std::size_t)i * n;
      for (std::size_t j = 0; j < i; ++j) {
        double v = col[j];                // == m[j + i*n] == m[i + j*n]
        if (ISNA(v) || ISNAN(v)) { v = na_value; ++local_imp; }
        int bl = dtt_format_value(buf, v, ndigits, fmt);
        row.append(buf, bl);
      }
      row.push_back('\n');
      (*imp)[i - base] = local_imp;
    }
  }
};

// [[Rcpp::export]]
double write_phylip_lower(NumericMatrix M, CharacterVector labels,
                          double na_value, std::string path,
                          int ndigits = 6) {
  const int n = M.nrow();
  if (M.ncol() != n)      stop("distance matrix must be square");
  if (labels.size() != n) stop("labels length must match matrix dimension");

  // 1. labels must be whitespace-free (PHYLIP parses the label as one token).
  //    Extract to a plain std::vector<std::string> here, in the R thread -- the
  //    parallel formatter must never touch the CharacterVector (R objects are
  //    not thread-safe to read).
  std::vector<std::string> labs;
  labs.reserve(n);
  for (int i = 0; i < n; ++i) {
    std::string lab = as<std::string>(labels[i]);
    if (lab.find_first_of(" \t\r\n") != std::string::npos)
      stop("label '%s' contains whitespace; PHYLIP labels may not", lab);
    labs.push_back(std::move(lab));
  }

  // 2. open the plain-text output (buffered stdio; no zlib involved).
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) stop("cannot open output file: %s", path);

  // 3. line 1: the taxon count.
  char hdr[32];
  int hl = std::snprintf(hdr, sizeof(hdr), "%d\n", n);
  std::fwrite(hdr, 1, (size_t)hl, f);

  const double* m = REAL(M);            // column-major: M(i,j) == m[i + j*n]

  // 4. per-row value format, e.g. " %.4f" (used only by the snprintf fallback).
  char fmt[16];
  std::snprintf(fmt, sizeof(fmt), " %%.%df", ndigits);

  // 5. format rows in parallel, one bounded chunk at a time, then fwrite each
  //    chunk in row order. Row i is independent (it only reads column i), so
  //    the output bytes are identical to a serial row-by-row writer.
  const int CHUNK = 4096;
  std::vector<std::string> bufs(std::min(CHUNK, n));
  std::vector<long long>   imp(std::min(CHUNK, n));
  double n_imputed = 0.0;

  for (int base = 0; base < n; base += CHUNK) {
    const int hi = std::min(base + CHUNK, n);
    RowFormatter worker(m, &labs, n, na_value, ndigits, fmt, base, &bufs, &imp);
    RcppParallel::parallelFor(base, hi, worker);
    for (int i = base; i < hi; ++i) {
      std::string& row = bufs[i - base];
      if (std::fwrite(row.data(), 1, row.size(), f) != row.size()) {
        std::fclose(f);
        stop("write failed while writing row %d", i + 1);
      }
      n_imputed += (double)imp[i - base];
    }
  }

  std::fclose(f);
  return n_imputed;
}
