// dtt_distance_scale.cpp
//
// Parallel kernel that reproduces, bit-for-bit, the cell-cell distance matrix
// of dtt_distance.R (dtt_distance_matrix).
//
// Reproduced quantity (per ordered pair of cells i, j):
//   shared = integrations recovered in BOTH cells (recovered = tape present and
//            not all-missing). If no shared integration -> NA.
//   For each shared integration t:
//       tape_distance = depth_i[t] + depth_j[t] - 2 * shared_edit_prefix
//   where depth = number of edited sites (id > 0) on the tape and
//   shared_edit_prefix = length of the leading run of sites that are the SAME
//   edit in both tapes (both id > 0 and equal), stopping at the first site that
//   is non-edit in either tape or a differing edit.
//   result = mean of tape_distance over shared integrations (a float64),
//   computed with R's two-pass mean.default algorithm so the value is
//   bit-identical to the reference: s = sum(d)/n; s += sum(d - s)/n.
//
// Encoding (built once in R, passed in here):
//   - code : IntegerVector length n*T*S, row-major [cell][tape][site], holding a
//            global edit-token id; id 0 == not-an-edit ("ETY"/"None"/"?"/etc.),
//            id > 0 == a distinct edit token. is_edit <=> id > 0.
//   - depth: IntegerVector length n*T, [cell][tape] = count of id>0 on that tape.
//   - mask : IntegerVector length n*W, the recovered-mask of each cell packed as
//            W = ceil(T/32) uint32 words (LSB-first); bit t set <=> tape t is
//            recovered in that cell.
//   Token equality is exact-string equality, captured by integer-id equality;
//   variable-width edit tokens (e.g. "AAG","GATG") are therefore safe.

#include <Rcpp.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;

struct DttWorker : public RcppParallel::Worker {
  const int* code;    // n*T*S
  const int* depth;   // n*T
  const unsigned int* mask;  // n*W
  const int n, T, S, W;
  double* out;        // n*n column-major (R matrix layout)

  DttWorker(const int* code_, const int* depth_, const unsigned int* mask_,
            int n_, int T_, int S_, int W_, double* out_)
    : code(code_), depth(depth_), mask(mask_), n(n_), T(T_), S(S_), W(W_),
      out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const unsigned int* mi = mask + i * W;
      const int* depth_i = depth + i * T;
      const int* code_i = code + i * (std::size_t)T * S;
      // Upper triangle (j >= i); mirror below. Each worker owns rows [begin,end)
      // and also writes the mirrored column entries for those rows -- but to
      // avoid cross-worker write races we write BOTH m[i,j] and m[j,i] only for
      // pairs where i is the row this worker owns, computing the full row j=0..n
      for (std::size_t j = 0; j < (std::size_t)n; ++j) {
        const unsigned int* mj = mask + j * W;
        int shared_count = 0;
        // R's mean.default does a two-pass mean over the per-tape distances:
        //   s = sum(d)/n;  s += sum(d - s)/n.
        // We reproduce that EXACTLY so the float64 result is bit-identical.
        // First pass: Sum of the per-tape distances (each a small integer,
        // exactly representable, so the sum is exact in double).
        double s1 = 0.0;
        const int* depth_j = depth + j * T;
        const int* code_j = code + j * (std::size_t)T * S;
        for (int w = 0; w < W; ++w) {
          unsigned int sh = mi[w] & mj[w];
          if (sh == 0u) continue;
          unsigned int bits = sh;
          while (bits) {
            int b = __builtin_ctz(bits);
            bits &= bits - 1;
            int t = w * 32 + b;
            ++shared_count;
            const int* ci = code_i + (std::size_t)t * S;
            const int* cj = code_j + (std::size_t)t * S;
            int prefix = 0;
            for (int s = 0; s < S; ++s) {
              int a = ci[s], c = cj[s];
              if (a > 0 && a == c) ++prefix; else break;
            }
            s1 += (double)(depth_i[t] + depth_j[t] - 2 * prefix);
          }
        }
        double val;
        if (shared_count == 0) {
          val = NA_REAL;
        } else {
          double nn = (double)shared_count;
          double s = s1 / nn;
          // Second pass: correction term, matching R's do_mean.
          double t = 0.0;
          for (int w = 0; w < W; ++w) {
            unsigned int sh = mi[w] & mj[w];
            if (sh == 0u) continue;
            unsigned int bits = sh;
            while (bits) {
              int b = __builtin_ctz(bits);
              bits &= bits - 1;
              int tt = w * 32 + b;
              const int* ci = code_i + (std::size_t)tt * S;
              const int* cj = code_j + (std::size_t)tt * S;
              int prefix = 0;
              for (int ss = 0; ss < S; ++ss) {
                int a = ci[ss], c = cj[ss];
                if (a > 0 && a == c) ++prefix; else break;
              }
              double d = (double)(depth_i[tt] + depth_j[tt] - 2 * prefix);
              t += d - s;
            }
          }
          val = s + t / nn;
        }
        // column-major store
        out[i + j * (std::size_t)n] = val;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix dtt_kernel(IntegerVector code, IntegerVector depth,
                         IntegerVector mask, int n, int T, int S, int W) {
  NumericMatrix out(n, n);
  DttWorker worker(code.begin(), depth.begin(),
                   (const unsigned int*)mask.begin(),
                   n, T, S, W, out.begin());
  RcppParallel::parallelFor(0, n, worker);
  return out;
}
