// nj_inmemory.cpp -- in-memory distance-matrix -> NJ tree handoff to decenttree,
// bypassing the text PHYLIP round-trip (write ~520 GB + serial re-parse ~7.8 h B1 /
// ~3.6 h B2). Reuses decenttree's OWN tested RapidNJ + parallel flat-array loader
// (headeronlydecenttree.cpp pattern), so we are not reimplementing NJ.
//
// Hosted in R via Rcpp: the distance matrix built by run_full_distance.R
// (`M`, a full n x n float64 R matrix, symmetric, diag 0, NA for unshared
// pairs) is handed straight to RapidNJ.
//
// Memory: R's M (8 n^2) + decenttree float32 working set (~12 n^2) coexist ->
//   B2 ~1.25 TB (fits ultramem-80 1.9 TB), B1 ~2.70 TB (needs 3.8 TB box).
//
// NA handling: decenttree cannot see NaN. R's NA_real_ IS a NaN bit pattern, so we
// impute NaN -> na_value IN PLACE on M's data (no R copy; caller discards M after),
// then call the tested loadDistancesFromFlatArray (parallel double->float cast).
//
// Precision: the production pipeline writes a PHYLIP at ndigits=4 (run_full_distance.R),
// so decenttree's CLI sees distances rounded to 4 decimals. To build the SAME tree as
// that pipeline (DTT distances are p/q with q<=11; full float64 vs 4-dp rounding flips
// many ties on real data), we round in place to `ndigits` decimals here too. ndigits<=0
// disables rounding (full precision -- used by the correctness oracle vs textbook NJ).
//
// CAUTION (in-place mutation): imputation + rounding overwrite M's own buffer, bypassing
// R's copy-on-write. Any R binding aliasing M (incl. one captured before this call) is
// mutated. Safe for the pipeline (M is discarded after), but never reuse M afterward.
//
// Threads: RapidNJ freezes threadCount = omp_get_max_threads() at construction, so we
// call omp_set_num_threads(nthreads) BEFORE constructing the object.

// mirror headeronlydecenttree.cpp's feature switches (keep it dependency-light)
#define USE_GZSTREAM            0
#define USE_PROGRESS_DISPLAY    0
// USE_VECTORCLASS_LIBRARY: deliberately left UNDEFINED, not defined to 0 --
// decenttree's upgma.h guards vectorclass's x86-only SIMD types with #ifdef
// (checks existence, not value), so `#define ... 0` would still pull it in
// and fail to compile on ARM64 (Apple Silicon has no SSE2). Every other use
// of this flag in decenttree's headers is a value check (#if), so leaving it
// undefined (treated as 0) disables vectorclass everywhere consistently.
#define DECENT_TREE             1

// [[Rcpp::plugins(cpp17)]]
// OpenMP flags are supplied explicitly via PKG_CXXFLAGS/PKG_LIBS by the R
// caller (see run_full_distance.R) rather than Rcpp's `plugins(openmp)`,
// which always emits bare `-fopenmp` -- correct for Linux gcc/clang, but
// unsupported by Apple clang (needs Homebrew libomp + different flags).
#include <Rcpp.h>
#include <cmath>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils/vectortypes.h" // StrVector (std::vector<std::string> subclass) for loadMatrix
#include "distancematrix.h"   // Matrix<T>::loadDistancesFromFlatArray
#include "rapidnj.h"          // StartTree::RapidNJ (== "NJ-R" == BoundingMatrix<NJFloat, NJMatrix<NJFloat>>)

using namespace Rcpp;

// [[Rcpp::export]]
void decenttree_nj_inmem(NumericMatrix M,
                         CharacterVector labels,
                         double         na_value,
                         int            nthreads,
                         int            precision,
                         std::string    out_path,
                         int            ndigits = 4) {
  const R_xlen_t n = labels.size();
  if (M.nrow() != n || M.ncol() != n) {
    stop("decenttree_nj_inmem: M must be n x n with n == length(labels).");
  }
  if (nthreads < 1) {
    stop("decenttree_nj_inmem: nthreads must be >= 1 (got %d).", nthreads);
  }
  // decenttree's loadMatrix ASSUMES 2 < names.size() (no runtime check of its
  // own), so guard here: NJ is undefined for < 3 taxa and a smaller n would let
  // decenttree read past a degenerate matrix. Fail with a catchable R error
  // instead of risking an abort/UB deep inside the library.
  if (n < 3) {
    stop("decenttree_nj_inmem: need at least 3 taxa (n == %d).", (int)n);
  }
  // The imputation value replaces every NA/NaN; if it is itself non-finite the
  // whole matrix would be poisoned and RapidNJ's minima become meaningless.
  if (!R_finite(na_value)) {
    stop("decenttree_nj_inmem: na_value must be finite (got non-finite value).");
  }

  // 1. thread count MUST be set before the algorithm object is constructed
  //    (RapidNJ snapshots omp_get_max_threads() in its ctor and never re-reads -nt).
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  // 2. impute NA/NaN -> na_value AND round to `ndigits` decimals, in place on R's own
  //    buffer (no copy). Rounding reproduces the pipeline's ndigits=4 PHYLIP so the tree
  //    matches the CLI on real p/q distances; ndigits<=0 leaves full precision.
  double* data = REAL(M);
  const R_xlen_t N = (R_xlen_t)n * (R_xlen_t)n;
  const bool   do_round = (ndigits > 0);
  const double scale    = do_round ? std::pow(10.0, ndigits) : 1.0;
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    double v = ISNAN(data[i]) ? na_value : data[i];
    if (do_round) v = std::nearbyint(v * scale) / scale;   // round-half-to-even, matches %.*f
    data[i] = v;
  }

  // 3. hand names + matrix to RapidNJ via the tested loadMatrix entry point.
  //    loadMatrix does the full, correct setup that a hand-rolled sequence misses:
  //    setSize -> clusters -> loadDistancesFromFlatArray (parallel double->float cast)
  //    -> calculateRowTotals (the per-row "U" sums NJ needs; omitting these yields a
  //    valid-looking but WRONG topology). R matrices are column-major, but M is
  //    symmetric with diag 0, so loadMatrix's row-major reading is identical to M.
  //    Guard labels: whitespace or a newick metacharacter would corrupt the tree
  //    (silent tip miscount), so reject them -- the CLI writer rejects whitespace too.
  StrVector names;
  names.reserve(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    std::string lab(labels[i]);
    if (lab.empty()) stop("decenttree_nj_inmem: tip label %d is empty.", (int)i + 1);
    if (lab.find_first_of(" \t\n\r(),:;") != std::string::npos) {
      stop("decenttree_nj_inmem: tip label '%s' contains whitespace or a newick "
           "metacharacter ( ) , : ; -- would corrupt the tree.", lab.c_str());
    }
    names.emplace_back(lab);
  }

  StartTree::RapidNJ r;
  r.loadMatrix(names, data);

  // 4. neighbour-join (multi-threaded clustering) and write the newick.
  r.constructTree();
  r.writeTreeFile(precision, out_path);
}
