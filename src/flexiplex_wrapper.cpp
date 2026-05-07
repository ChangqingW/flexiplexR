#include "main-functions/flexiplex.h"
#include <Rcpp.h>

//' Rcpp port of flexiplex
//'
//' @description demultiplex reads with flexiplex. For a detailed description
//' see documentation for the original flexiplex:
//' https://davidsongroup.github.io/flexiplex
//'
//' @param r_segments R list of FlexiplexSegment S4 objects
//' @param r_barcode_groups R list of FlexiplexGroup S4 objects
//' @param max_flank_editdistance maximum edit distance for flanking sequences
//' @param reads_in input FASTQ or FASTA file path(s)
//' @param reads_out output file for demultiplexed reads
//' @param stats_out output file for per-read demultiplex statistics
//' @param bc_out output file for barcode counts
//' @param reverseCompliment whether to reverse complement reads after demultiplexing
//' @param n_threads number of threads
//' @return IntegerVector with read counts (total, demultiplexed, single match, chimera)
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector flexiplex(
  Rcpp::List r_segments,
  Rcpp::List r_barcode_groups,
  int max_flank_editdistance,
  Rcpp::StringVector reads_in,
  Rcpp::String reads_out,
  Rcpp::String stats_out,
  Rcpp::String bc_out, bool reverseCompliment, int n_threads) {

  return flexiplex_cpp(
    r_segments,
    r_barcode_groups,
    max_flank_editdistance,
    reads_in,
    reads_out,
    stats_out,
    bc_out, reverseCompliment, n_threads
  );
}
