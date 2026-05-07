#ifndef FLEXIPLEX_H
#define FLEXIPLEX_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector flexiplex_cpp(
  Rcpp::List r_segments,
  Rcpp::List r_barcode_groups,
  int max_flank_editdistance,
  Rcpp::StringVector reads_in,
  Rcpp::String reads_out,
  Rcpp::String stats_out,
  Rcpp::String bc_out, bool reverseCompliment, int n_threads);

// Note: `edit_distance()` is implemented as a `static inline` helper inside
// `flexiplex.cpp` and is intentionally not declared in this header to avoid
// ODR/linkage conflicts.

#endif // FLEXIPLEX_H
