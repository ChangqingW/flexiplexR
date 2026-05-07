#include "flexiplex.h"
// Copyright 2022 Nadia Davidson
// This program is distributed under the MIT License.
// We also ask that you cite this software in publications
// where you made use of it for any part of the data analysis.

#include <Rcpp.h>
#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <istream>
#include <numeric>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../utility/edlib-1.2.7/edlib.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "zlib.h"
// [[Rcpp::plugins(cpp17)]]

#ifndef GZKSEQ
#define GZKSEQ
KSEQ_INIT(gzFile, gzread)
#endif

// Append .1 to version for dev code, remove for release
// e.g. 1.00.1 (dev) goes to 1.01 (release)
const static std::string VERSION = "1.02.6";

enum SegmentType { FIXED, MATCHED, MATCHED_SPLIT, RANDOM };

struct Segment {
    SegmentType type;
    std::string pattern;
    std::string name; // e.g., "BC1", "UMI", "Primer1"
    std::string bc_list_name; // Only for MATCHED, refers to a key in known_barcodes_map
    int buffer_size;            // for MATCHED, MATCHED_SPLIT
    int max_edit_distance;      // for MATCHED, MATCHED_SPLIT
};

struct BarcodeGroup {
    std::string name;
    int max_edit_distance;
    std::vector<size_t> segment_indices;
};

// complement nucleotides - used to reverse complement string
char complement(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//Inplace reverse complement
void reverse_complement(std::string & seq){
  std::reverse(seq.begin(),seq.end());
  std::transform(seq.begin(),seq.end(),seq.begin(),complement);
}

std::string reverse_complement_copy(const std::string & seq){
  std::string rc_seq = seq;
  std::reverse(rc_seq.begin(),rc_seq.end());
  std::transform(rc_seq.begin(),rc_seq.end(),rc_seq.begin(),complement);
  return rc_seq;
}

// Holds the found barcode and associated information
struct Barcode {
  std::map<std::string, std::string> features; // Map of segment name to extracted sequence
  int flank_editd = 100;
  int flank_start = -1;
  int flank_end = -1; // inclusive
  bool found_all_matched_segments = true;
};

struct SearchResult {
  std::string read_id;
  std::string qual_scores;
  std::string line;
  std::string rev_line;
  std::vector<Barcode> vec_bc_for;
  std::vector<Barcode> vec_bc_rev;
  int count = 0;
  bool chimeric = false;
};

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumed to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd) {

  const std::size_t len1 = s1.size() + 1, len2 = s2.size() + 1;
  const char *s1_c = s1.c_str(); const char *s2_c = s2.c_str();

  // Reuse DP buffer per-thread to reduce allocator churn.
  thread_local std::vector<unsigned int> dist_holder;
  dist_holder.assign(len1 * len2, 0u);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0]=0; //[0][0]
  for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
  for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];

  int best = len2;
  end = len1 - 1;

  //loop over the distance matrix elements and calculate running distance
  for (unsigned int j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (unsigned int i = 1; i < len1; ++i) {
      const int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1; // are the bases the same?

      // if yes, no need to increment distance
      if (sub == 0) {
        dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
      } else {
        dist_holder[i * len2 + j] =
            std::min({dist_holder[(i - 1) * len2 + j] + 1,
                      dist_holder[i * len2 + (j - 1)] + 1,
                      dist_holder[(i - 1) * len2 + (j - 1)] + 1});
      }

      if (dist_holder[i * len2 + j] <= (unsigned)max_editd) {
        any_below_threshold = true;
      }

      // if this is the last row in j
      if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
        // check if this is the best running score
        best = dist_holder[i * len2 + j];
        end = i; // update the end position of alignment
      }
    }

    // early exit to save time.
    if (!any_below_threshold) {
      return 100;
    }
  }

  return best; // return edit distance
}

// Forward declarations
bool align_read_to_pattern(const std::string &seq,
                           const std::vector<Segment> &segments,
                           int global_flank_max_editd,
                           Barcode &barcode_result,
                           std::vector<int> &read_to_segment_starts);
void refine_matched_segments(const std::string &seq,
                             const std::vector<Segment> &segments,
                             const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                             const std::vector<int> &read_to_segment_starts,
                             Barcode &barcode_result,
                             std::vector<int> &refined_segment_starts,
                             std::vector<int> &refined_segment_ends);
void refine_split_segments(const std::string &seq,
                           const std::vector<Segment> &segments,
                           const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                           const std::unordered_map<std::string, BarcodeGroup> *group_map,
                           const std::vector<int> &read_to_segment_starts,
                           Barcode &barcode_result,
                           std::vector<int> &refined_segment_starts,
                           std::vector<int> &refined_segment_ends);
void extract_random_segments(const std::string &seq,
                             const std::vector<Segment> &segments,
                             const std::vector<int> &read_to_segment_starts,
                             const std::vector<int> &refined_segment_starts,
                             const std::vector<int> &refined_segment_ends,
                             Barcode &barcode_result);

// Given a string 'seq' search for substring with primer and polyT sequence followed by
// a targeted search in the region for barcode
// Sequence seearch is performed using edlib

Barcode extract_features(std::string & seq,
                         const std::vector<Segment> &segments,
                         const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                         const std::unordered_map<std::string, BarcodeGroup> *group_map,
                         int global_flank_max_editd) {

  Barcode barcode_result;
  barcode_result.found_all_matched_segments = true;
  barcode_result.flank_editd = 100;

  std::vector<int> read_to_segment_starts;

  // 1. Align reads to get approximate positions
  bool aligned = align_read_to_pattern(seq, segments, global_flank_max_editd, barcode_result, read_to_segment_starts);

  if (!aligned) return barcode_result;

  // Storage for refined positions from matched segments
  std::vector<int> refined_segment_starts(segments.size(), -1);
  std::vector<int> refined_segment_ends(segments.size(), -1);

  // 2. Process MATCHED segments (Single barcodes)
  refine_matched_segments(seq, segments, known_barcodes_map, read_to_segment_starts, barcode_result, refined_segment_starts, refined_segment_ends);

  // 3. Process MATCHED_SPLIT segments (Grouped/Split barcodes)
  if (group_map && !group_map->empty()) {
    refine_split_segments(seq, segments, known_barcodes_map, group_map, read_to_segment_starts, barcode_result, refined_segment_starts, refined_segment_ends);
  }

  if (!barcode_result.found_all_matched_segments) {
    return barcode_result;
  }

  // 4. Process RANDOM segments using refined anchors
  extract_random_segments(seq, segments, read_to_segment_starts, refined_segment_starts, refined_segment_ends, barcode_result);

  return barcode_result;
}

// Helper: Align read to the concatenated pattern using Edlib
bool align_read_to_pattern(const std::string &seq,
                           const std::vector<Segment> &segments,
                           int global_flank_max_editd,
                           Barcode &barcode_result,
                           std::vector<int> &read_to_segment_starts) {

  // Edlib Configuration 
  EdlibEqualityPair additionalEqualities[32] = {
    {'R', 'A'}, {'R', 'G'}, {'K', 'G'}, {'K', 'T'},
    {'S', 'G'}, {'S', 'C'}, {'Y', 'C'}, {'Y', 'T'},
    {'M', 'A'}, {'M', 'C'}, {'W', 'A'}, {'W', 'T'},
    {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
    {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
    {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'},
    {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
  };
  EdlibAlignConfig edlibConf = {global_flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 32};

  // 1. Construct Search Template
  // Optimization Note: This construction happens every call. If segments are static, 
  // the template string and lengths could be pre-calculated.
  std::string search_string_template;
  std::vector<long unsigned int> segment_lengths;

  search_string_template.reserve(256);
  segment_lengths.reserve(segments.size());

  for (const auto &segment : segments) {
    search_string_template += segment.pattern;
    segment_lengths.push_back(segment.pattern.length());
  }

  if (seq.length() < search_string_template.length()) return false;

  // 2. Perform Alignment
  EdlibAlignResult result = edlibAlign(search_string_template.c_str(), search_string_template.length(), seq.c_str(), seq.length(), edlibConf);

  if (result.status != EDLIB_STATUS_OK || result.numLocations == 0) {
    edlibFreeAlignResult(result);
    return false;
  }

  barcode_result.flank_editd = result.editDistance;
  barcode_result.flank_start = result.startLocations[0];
  barcode_result.flank_end = result.endLocations[0];

  // 3. Map alignment to segment positions
  std::vector<long unsigned int> segment_template_ends;
  segment_template_ends.resize(segment_lengths.size());
  std::partial_sum(segment_lengths.begin(), segment_lengths.end(), segment_template_ends.begin());

  read_to_segment_starts.reserve(segment_template_ends.size() + 1);
  read_to_segment_starts.emplace_back(barcode_result.flank_start);

  int i_read = barcode_result.flank_start;
  int i_pattern = 0;
  size_t i_segment = 0;

  // Avoid copying the alignment to a temporary vector.
  for (int ai = 0; ai < result.alignmentLength; ++ai) {
    const unsigned char value = result.alignment[ai];
    if (value != EDLIB_EDOP_INSERT) i_read++;
    if (value != EDLIB_EDOP_DELETE) i_pattern++;

    if (i_segment < segment_template_ends.size() &&
        i_pattern >= (int)segment_template_ends[i_segment]) {
      read_to_segment_starts.emplace_back(i_read);
      i_segment++;
    }
  }

  edlibFreeAlignResult(result);
  return true;
}

// Helper: Process direct MATCHED segments (Pass 1)
void refine_matched_segments(const std::string &seq,
                             const std::vector<Segment> &segments,
                             const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                             const std::vector<int> &read_to_segment_starts,
                             Barcode &barcode_result,
                             std::vector<int> &refined_segment_starts,
                             std::vector<int> &refined_segment_ends) {

  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment& s = segments[i];
    if (s.type != MATCHED) continue;

    const int segment_read_start = read_to_segment_starts[i];
    const int segment_read_end =
      (i + 1 < read_to_segment_starts.size())
      ? read_to_segment_starts[i+1]
      : barcode_result.flank_end;

    const std::unordered_set<std::string>* current_bclist = nullptr;
    if (!s.bc_list_name.empty()) {
      auto it = known_barcodes_map->find(s.bc_list_name);
      if (it != known_barcodes_map->end()) current_bclist = &(it->second);
    } else if (known_barcodes_map->count("global")) {
      current_bclist = &(known_barcodes_map->at("global"));
    }

    if (current_bclist && !current_bclist->empty()) {
      const int search_start = std::max(0, segment_read_start - s.buffer_size);
      const int search_end = std::min((int)seq.length(), segment_read_end + s.buffer_size);
      const std::string search_region = seq.substr(search_start, search_end - search_start);

      unsigned int best_edit_distance = 100;
      unsigned int best_end_distance = 0;
      std::string best_barcode_match = "";
      bool current_segment_unambiguous = false;

      for (const auto& known_bc : *current_bclist) {
        unsigned int editDistance, endDistance;
        editDistance = edit_distance(search_region, known_bc, endDistance, s.max_edit_distance);

        if (editDistance == best_edit_distance) {
          current_segment_unambiguous = false;
        } else if (editDistance < best_edit_distance && editDistance <= s.max_edit_distance) {
          current_segment_unambiguous = true;
          best_edit_distance = editDistance;
          best_barcode_match = known_bc;
          best_end_distance = endDistance;
          if (editDistance == 0) break; 
        }
      }
      if (best_edit_distance <= s.max_edit_distance && current_segment_unambiguous) {
        barcode_result.features[s.name] = best_barcode_match;
        refined_segment_ends[i] = search_start + best_end_distance - 1;
        refined_segment_starts[i] = refined_segment_ends[i] - best_barcode_match.length() + 1;
      } else {
        barcode_result.found_all_matched_segments = false;
      }
    } else {
      // Discovery Mode (Fallback when no list provided)
      const int len = segment_read_end - segment_read_start;
      if (len > 0 && segment_read_start + len <= (int)seq.length()) {
        barcode_result.features[s.name] = seq.substr(segment_read_start, len);
        refined_segment_starts[i] = segment_read_start;
        refined_segment_ends[i] = segment_read_end - 1;
      } else {
        barcode_result.features[s.name] = "";
        refined_segment_starts[i] = segment_read_start;
        refined_segment_ends[i] = segment_read_start - 1;
      }
    }
  }
}

// Helper: Process MATCHED_SPLIT segments (Pass 2)
void refine_split_segments(const std::string &seq,
                           const std::vector<Segment> &segments,
                           const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                           const std::unordered_map<std::string, BarcodeGroup> *group_map,
                           const std::vector<int> &read_to_segment_starts,
                           Barcode &barcode_result,
                           std::vector<int> &refined_segment_starts,
                           std::vector<int> &refined_segment_ends) {

  for (size_t i = 0; i < segments.size(); ++i) {
    if (segments[i].type != MATCHED_SPLIT) continue;
    if (refined_segment_starts[i] != -1) continue; // Already processed as part of a group

    const std::string group_name = segments[i].bc_list_name;
    if (group_map->find(group_name) == group_map->end()) {
      Rcpp::stop("Error: Undefined group " + group_name + ".\n");
    }

    const BarcodeGroup& bg = group_map->at(group_name);
    const std::vector<size_t>& split_group_indices = bg.segment_indices;

    std::vector<std::string> part_seqs;
    std::vector<int> part_starts;
    part_seqs.reserve(split_group_indices.size());
    part_starts.reserve(split_group_indices.size());
    bool possible_to_extract = true;

    for (size_t idx : split_group_indices) {
      const Segment &s = segments[idx];
      const int segment_read_start = read_to_segment_starts[idx];
      const int segment_read_end =
        (idx + 1 < read_to_segment_starts.size())
        ? read_to_segment_starts[idx + 1]
        : barcode_result.flank_end;

      const int search_start = std::max(0, segment_read_start - s.buffer_size);
      const int search_end = std::min((int)seq.length(), segment_read_end + s.buffer_size);

      if (search_end <= search_start) {
        possible_to_extract = false;
        break;
      }

      part_seqs.emplace_back(seq.substr(search_start, search_end - search_start));
      part_starts.emplace_back(search_start);
    }

    if (!possible_to_extract) {
      barcode_result.found_all_matched_segments = false;
      continue;
    }

    // Lookup group barcode list. If not present or empty, fall back to discovery mode.
    const std::unordered_set<std::string>* current_bclist = nullptr;
    auto it = known_barcodes_map->find(group_name);
    if (it != known_barcodes_map->end()) current_bclist = &(it->second);

    if (current_bclist && !current_bclist->empty()) {
      unsigned int best_edit_distance = 100;
      std::string best_barcode_match = "";
      std::vector<int> best_part_starts_in_read(split_group_indices.size(), 0);
      bool current_group_unambiguous = false;

      std::vector<size_t> known_split_offsets;
      known_split_offsets.reserve(split_group_indices.size());
      size_t running_offset = 0;
      for (size_t idx : split_group_indices) {
        known_split_offsets.push_back(running_offset);
        running_offset += segments[idx].pattern.length();
      }

      for (const auto& known_bc : *current_bclist) {
        unsigned int total_edit_distance = 0;
        std::vector<int> current_part_starts;
        current_part_starts.reserve(split_group_indices.size());
        bool possible_match = true;

        for (size_t k = 0; k < split_group_indices.size(); ++k) {
          const size_t idx = split_group_indices[k];
          const int offset = known_split_offsets[k];
          int len = segments[idx].pattern.length();

          if (offset + len > (int)known_bc.length()) len = (int)known_bc.length() - offset;
          if (len <= 0) {
            possible_match = false;
            break;
          }

          const std::string known_part = known_bc.substr(offset, len);
          unsigned int endDist = 0;
          const unsigned int part_ed =
            edit_distance(part_seqs[k], known_part, endDist,
                          segments[idx].max_edit_distance);

          if (part_ed > (unsigned)segments[idx].max_edit_distance) {
            possible_match = false;
            break;
          }
          total_edit_distance += part_ed;
          current_part_starts.push_back(part_starts[k] + endDist - known_part.length());
        }

        if (!possible_match)
          continue;

        if (total_edit_distance == best_edit_distance) {
          current_group_unambiguous = false;
        } else if (total_edit_distance < best_edit_distance) {
          current_group_unambiguous = true;
          best_edit_distance = total_edit_distance;
          best_barcode_match = known_bc;
          best_part_starts_in_read = current_part_starts;
        }
      }

      if (best_edit_distance <= (unsigned)bg.max_edit_distance && current_group_unambiguous) {
        size_t running_offset = 0;
        for (size_t k = 0; k < split_group_indices.size(); ++k) {
          const size_t idx = split_group_indices[k];
          int len = segments[idx].pattern.length();
          if (running_offset + len > best_barcode_match.length())
            len = best_barcode_match.length() - running_offset;

          barcode_result.features[segments[idx].name] =
            best_barcode_match.substr(running_offset, len);
          refined_segment_starts[idx] = best_part_starts_in_read[k];
          refined_segment_ends[idx] = best_part_starts_in_read[k] + len - 1;
          running_offset += segments[idx].pattern.length();
        }
        barcode_result.features[group_name] = best_barcode_match;
      } else {
        barcode_result.found_all_matched_segments = false;
      }

    } else {
      // Discovery mode for MATCHED_SPLIT:
      // No known list for the group, so extract each part based on the approximate
      // alignment-derived boundaries and concatenate to a group-level feature.
      std::string group_concat;
      group_concat.reserve(256);

      for (size_t k = 0; k < split_group_indices.size(); ++k) {
        const size_t idx = split_group_indices[k];
        const Segment& s = segments[idx];

        int segment_read_start = read_to_segment_starts[idx];
        int segment_read_end =
          (idx + 1 < read_to_segment_starts.size())
          ? read_to_segment_starts[idx + 1]
          : barcode_result.flank_end;

        int extract_start = std::max(0, segment_read_start);
        int extract_end = std::min((int)seq.length(), segment_read_end);

        // Anchor left boundary to previous refined split-part if we have it.
        if (k > 0) {
          const size_t prev_idx = split_group_indices[k - 1];
          if (refined_segment_ends[prev_idx] != -1) {
            extract_start = std::max(extract_start, refined_segment_ends[prev_idx] + 1);
          }
        }
        // Anchor right boundary to next refined split-part if we have it.
        if (k + 1 < split_group_indices.size()) {
          const size_t next_idx = split_group_indices[k + 1];
          if (refined_segment_starts[next_idx] != -1) {
            extract_end = std::min(extract_end, refined_segment_starts[next_idx]);
          }
        }

        // If still unresolved, fall back to expected pattern length, starting at the approximate position.
        if (extract_end <= extract_start) {
          extract_start = std::max(0, segment_read_start);
          extract_end = std::min((int)seq.length(), extract_start + (int)s.pattern.length());
        }

        if (extract_end < extract_start) extract_end = extract_start;

        const std::string part = seq.substr(extract_start, extract_end - extract_start);
        barcode_result.features[s.name] = part;
        refined_segment_starts[idx] = extract_start;
        refined_segment_ends[idx] = extract_start + (int)part.length() - 1;

        group_concat += part;
      }

      barcode_result.features[group_name] = group_concat;
    }
  }
}

// Helper: Process RANDOM segments (Pass 3)
void extract_random_segments(const std::string &seq,
                             const std::vector<Segment> &segments,
                             const std::vector<int> &read_to_segment_starts,
                             const std::vector<int> &refined_segment_starts,
                             const std::vector<int> &refined_segment_ends,
                             Barcode &barcode_result) {

  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment &s = segments[i];
    if (s.type == MATCHED || s.type == MATCHED_SPLIT) continue;

    int extract_start = read_to_segment_starts[i];
    int extract_end = (i + 1 < read_to_segment_starts.size())
      ? read_to_segment_starts[i + 1]
      : barcode_result.flank_end;

    if (s.type == RANDOM) {
      // Anchor to previous matched segment if available
      if (i > 0 &&
        (segments[i - 1].type == MATCHED || segments[i - 1].type == MATCHED_SPLIT) &&
          refined_segment_ends[i - 1] != -1) {
        extract_start = refined_segment_ends[i - 1] + 1;
        extract_end = extract_start + (int)s.pattern.length();
      }
      // Anchor to next matched segment if available
      else if (i + 1 < segments.size() &&
        (segments[i + 1].type == MATCHED || segments[i + 1].type == MATCHED_SPLIT) &&
          refined_segment_starts[i + 1] != -1) {
        extract_end = refined_segment_starts[i + 1];
        extract_start = extract_end - (int)s.pattern.length();
      }

      int actual_extract_start = std::max(0, extract_start);
      int actual_extract_end = std::min((int)seq.length(), extract_end);

      if (actual_extract_end < actual_extract_start) actual_extract_end = actual_extract_start;

      barcode_result.features[s.name] =
        seq.substr(actual_extract_start, actual_extract_end - actual_extract_start);
    }
  }
}

std::vector<Barcode> big_barcode_search(
  const std::string &sequence,
  const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
  const std::unordered_map<std::string, BarcodeGroup> *group_map,
  int global_flank_max_editd,
  const std::vector<Segment> &segments) {

  std::vector<Barcode> return_vec;
  std::string masked_sequence = sequence; // Work on a copye

  while (true) {
    Barcode result = extract_features(masked_sequence, segments, known_barcodes_map,
                                      group_map, global_flank_max_editd);

    if (result.flank_editd <= global_flank_max_editd && result.found_all_matched_segments) {
      return_vec.emplace_back(result);

      // Mask the found region to prevent re-finding it
      // result.flank_end is inclusive, so length is end - start + 1
      const int match_length = result.flank_end - result.flank_start + 1;

      if (match_length > 0) {
        masked_sequence.replace(result.flank_start, match_length, std::string(match_length, 'X'));
      } else {
        break; // Should not happen for valid match, but prevents infinite loop
      }
    } else {
      break;
    }
  }
  return return_vec;
}

static void write_bgzfstring(BGZF *bgzf, const std::string &str) {
  const char *data = str.c_str();
  size_t len = str.length();

  int bytes_written = bgzf_write(bgzf, data, len);
  if (bytes_written != (int)str.size()) {
    Rcpp::stop("BGZF write error while writing string: %s\n, expected %zu bytes, wrote %d bytes",
              str.c_str(), str.size(), bytes_written);
  }
}

// print information about barcodes
void print_stats(const std::string &read_id, const std::vector<Barcode> &vec_bc,
                        BGZF *bgzf, const std::vector<Segment> &segments) {
  for (const auto &bc : vec_bc) {
    std::string line = read_id;
    for (const auto &s : segments) {
      if (s.type == FIXED)
        continue;
      auto it = bc.features.find(s.name);
      line += "\t" + (it != bc.features.end() ? it->second : std::string());
    }
    line += "\t" + std::to_string(bc.flank_editd) + "\n";
    write_bgzfstring(bgzf, line);
  }
}

void print_line(const std::string &id, const std::string &read,
                       const std::string &quals, bool reverse_complement_out,
                       BGZF *bgzf) {

  const char delimiter = quals.empty() ? '>' : '@';

  if (reverse_complement_out) {
    std::string rev_seq_lines =
        delimiter + id + "\n" + reverse_complement_copy(read) + "\n";
    write_bgzfstring(bgzf, rev_seq_lines);
    if (!quals.empty()) {
      std::string rev_qual_lines =
          "+\n" + std::string(quals.rbegin(), quals.rend()) + "\n";
      write_bgzfstring(bgzf, rev_qual_lines);
    }
  } else {
    std::string seq_lines = delimiter + id + "\n" + read + "\n";
    write_bgzfstring(bgzf, seq_lines);
    if (!quals.empty()) {
      std::string qual_lines = "+\n" + quals + "\n";
      write_bgzfstring(bgzf, qual_lines);
    }
  }
}

std::string compose_new_id(
    const std::string &read_id, const Barcode &bc, int which, int total,
    bool chimeric, const std::vector<Segment> &segments,
    const std::unordered_map<std::string, BarcodeGroup> &group_map) {

  std::ostringstream ss_suffix;
  ss_suffix << which << "of" << total;
  if (chimeric)
    ss_suffix << "_C";

  std::ostringstream id;

  // print the grouped barcodes as a whole first
  std::string delim = "";
  for (const auto &g : group_map) {
    auto it = bc.features.find(g.first);
    if (it != bc.features.end()) {
      id << delim << it->second;
    } else {
      id << delim << "NA";
    }
    delim = "-";
  }

  // print barcodes
  for (const auto &s : segments) {
    if (s.type != MATCHED)
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      id << delim << it->second;
    } else {
      id << delim << "NA";
    }
    delim = "-";
  }

  // Add UB (if present, concatonate all UB segments if multiple present)
  std::string ub = "";
  for (const auto &s : segments) {
    if (s.type != RANDOM)
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      ub += it->second;
    }
  }
  if (!ub.empty()) {
    id << "_" << ub;
  } else {
    id << "_" << "NA";
  }

  // Add CB tag (assume always present)
  id << "#" << read_id << ss_suffix.str() << "\tCB:Z:";
  delim = "";
  for (const auto &g : group_map) {
    auto it = bc.features.find(g.first);
    if (it != bc.features.end()) {
      id << delim << it->second;
      delim = ",";
    }
  }
  for (const auto &s : segments) {
    if (s.type != MATCHED)
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      id << delim << it->second;
      delim = ",";
    }
  }

  // add UB tag if exist
  if (!ub.empty()) {
    id << "\tUB:Z:" << ub;
  }

  return id.str();
}

void print_read(const std::string &read_id, const std::string &read,
                       const std::string &qual,
                       const std::vector<Barcode> &vec_bc, BGZF *bgzf,
                       bool trim_barcodes, bool chimeric,
                       bool reverse_complement_out,
                       const std::vector<Segment> &segments,
                       const std::unordered_map<std::string, BarcodeGroup>&group_map) {

  const size_t vec_size = vec_bc.size();

  for (int b = 0; b < vec_size; b++) {
    const Barcode &bc = vec_bc.at(b);

    if (bc.flank_end < 0) {
      continue;
    }

    std::string new_read_id =
      compose_new_id(read_id, bc, b + 1, vec_size, chimeric, segments, group_map);

    int read_start = bc.flank_end + 1;
    // work out the start and end base in case multiple barcodes
    int next_barcode_start = read.length();
    for (int k = 0; k < vec_size; k++) {
      if (b == k) continue;
      if (vec_bc.at(k).flank_start >= read_start) {
        next_barcode_start = std::min(next_barcode_start, vec_bc.at(k).flank_start);
      }
    }
    int read_length = next_barcode_start - read_start;
    read_length = std::max(0, read_length);

    std::string qual_new;
    if (!qual.empty()) {
      if (read_start + read_length > (int)qual.length()) {
        Rcpp::warning("sequence and quality lengths diff for read:  %s. Ignoring read.\n",
                     read_id.c_str());
        return;
      }
      qual_new = qual.substr(read_start, read_length);
    }
    std::string read_new = read.substr(read_start, read_length);

    if (b == 0 && !trim_barcodes) {
      new_read_id = read_id;
      read_new = read;
      qual_new = qual;
    }

    if (read_new.empty()) {
      continue;
    }

    print_line(new_read_id, read_new, qual_new, reverse_complement_out, bgzf);
  }
}

static void search_read(
    std::vector<SearchResult> &reads,
    const std::unordered_map<std::string, std::unordered_set<std::string>>
        &known_barcodes_map,
    const std::unordered_map<std::string, BarcodeGroup> &group_map,
    int flank_edit_distance,
    const std::vector<Segment> &segments) {

  for (auto &read : reads) {
    read.vec_bc_for = big_barcode_search(read.line, &known_barcodes_map,
                                        &group_map, flank_edit_distance,
                                        segments);

    read.rev_line = reverse_complement_copy(read.line);

    read.vec_bc_rev = big_barcode_search(read.rev_line, &known_barcodes_map,
                                        &group_map, flank_edit_distance,
                                        segments);

    read.count = read.vec_bc_for.size() + read.vec_bc_rev.size();
    read.chimeric = !read.vec_bc_for.empty() && !read.vec_bc_rev.empty();
  }
}

static bool file_exists(const std::string &filename) {
  std::ifstream infile(filename);
  return infile.good();
}

std::unordered_set<std::string> load_barcode_list(const std::string &filename) {
  std::unordered_set<std::string> bclist;
  std::ifstream bc_file(filename);
  if (bc_file.good()) {
    Rcpp::Rcout << "Loading known barcodes from " << filename << "\n";
    std::string line;
    std::string bc;
    while (getline(bc_file, line)) {
      std::istringstream line_stream(line);
      line_stream >> bc; // works for whitespace / tab delimited files
      if (!bc.empty())
        bclist.insert(bc);
    }
    Rcpp::Rcout << "Number of known barcodes: " << bclist.size() << "\n";
  } else {
    Rcpp::stop("Error: Unable to open barcode list file: " + filename + "\n");
  }
  return bclist;
}

// Helpers for R interface
SegmentType segment_type_from_string(const std::string &type_str) {
  if (type_str == "FIXED") {
    return FIXED;
  } else if (type_str == "MATCHED") {
    return MATCHED;
  } else if (type_str == "MATCHED_SPLIT") {
    return MATCHED_SPLIT;
  } else if (type_str == "RANDOM") {
    return RANDOM;
  } else {
    Rcpp::stop("Unknown segment type: " + type_str);
  }
  return FIXED; // Unreachable, but silences compiler warning
}

Segment s4_to_segment(const Rcpp::S4 &obj) {

  Segment seg;
  seg.type = segment_type_from_string(Rcpp::as<std::string>(obj.slot("type")));
  seg.pattern = Rcpp::as<std::string>(obj.slot("pattern"));
  seg.name = Rcpp::as<std::string>(obj.slot("name"));

  Rcpp::CharacterVector r_bc_list_name = obj.slot("bc_list_name");
  if (r_bc_list_name.size() > 0 && !Rcpp::CharacterVector::is_na(r_bc_list_name[0])) {
    seg.bc_list_name = Rcpp::as<std::string>(r_bc_list_name[0]);
  }

  Rcpp::IntegerVector r_buffer_size = obj.slot("buffer_size");
  if (r_buffer_size.size() > 0 && !Rcpp::IntegerVector::is_na(r_buffer_size[0])) {
    seg.buffer_size = r_buffer_size[0];
  }

  Rcpp::IntegerVector r_max_edit_distance = obj.slot("max_edit_distance");
  if (r_max_edit_distance.size() > 0 && !Rcpp::IntegerVector::is_na(r_max_edit_distance[0])) {
    seg.max_edit_distance = r_max_edit_distance[0];
  }

  return seg;
}

BarcodeGroup s4_to_barcode_group(const Rcpp::S4 &obj) {
  BarcodeGroup bg;
  bg.name = Rcpp::as<std::string>(obj.slot("name"));
  bg.max_edit_distance = Rcpp::as<int>(obj.slot("max_edit_distance"));
  return bg;
}

//' Rcpp port of flexiplex
//'
//' @description demultiplex reads with flexiplex, for detailed description, see
//' documentation for the original flexiplex:
//' https://davidsongroup.github.io/flexiplex
//'
//' @param reads_in Input FASTQ or FASTA file
//' @param barcodes_file barcode allow-list file
//' @param bc_as_readid bool, whether to add the demultiplexed barcode to the
//' read ID field
//''
//' @param r_segments R list of Segment S4 objects
//' @param r_barcode_groups R list of BarcodeGroup S4 objects
//' @param reads_out output file for demultiplexed reads
//' @param stats_out output file for demultiplexed stats
//' @param n_threads number of threads to be used during demultiplexing
//' @param reverseCompliment bool, whether to reverse complement the reads after demultiplexing
//' @param bc_out WIP
//' @return integer return value. 0 represents normal return.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector flexiplex_cpp(
  Rcpp::List r_segments,
  Rcpp::List r_barcode_groups,
  int max_flank_editdistance,
  Rcpp::StringVector reads_in,
  Rcpp::String reads_out,
  Rcpp::String stats_out,
  Rcpp::String bc_out, bool reverseCompliment, int n_threads) {

  std::ios_base::sync_with_stdio(false);

  // flexiplex -i option
  const bool remove_barcodes = true;

  // Parse arguments from R
  std::unordered_map<std::string, std::unordered_set<std::string>>
      known_barcodes_map;
  // segments
  std::vector<Segment> segments;
  segments.reserve(r_segments.size());
  for (const auto &r_seg : r_segments) {
    Segment seg = s4_to_segment(Rcpp::as<Rcpp::S4>(r_seg));
    segments.push_back(seg);
    if (seg.type == MATCHED && !seg.bc_list_name.empty() &&
      known_barcodes_map.find(seg.bc_list_name) == known_barcodes_map.end()) {
      known_barcodes_map[seg.bc_list_name] = load_barcode_list(seg.bc_list_name);
    }
  }
  // groups
  std::unordered_map<std::string, BarcodeGroup> group_map;
  for (const auto &r_grp : r_barcode_groups) {
    BarcodeGroup bg = s4_to_barcode_group(Rcpp::as<Rcpp::S4>(r_grp));
    group_map[bg.name] = bg;
    Rcpp::CharacterVector r_bc_list_file = Rcpp::as<Rcpp::S4>(r_grp).slot("bc_list_name");
    if (r_bc_list_file.size() > 0 && !Rcpp::CharacterVector::is_na(r_bc_list_file[0])) {
      known_barcodes_map[bg.name] = load_barcode_list(Rcpp::as<std::string>(r_bc_list_file[0]));
    } else {
      Rcpp::stop("Error: Barcode group " + bg.name + " requires a barcode list file.\n");
    }
  }
  // populate group segment indices
  for (size_t i = 0; i < segments.size(); ++i) {
    if (segments[i].type == MATCHED_SPLIT) {
      if (group_map.find(segments[i].bc_list_name) != group_map.end()) {
        group_map[segments[i].bc_list_name].segment_indices.push_back(i);
      } else {
        Rcpp::stop("Error: Undefined group " + segments[i].bc_list_name + ".\n");
      }
    }
  }

  Rcpp::Rcout << "FLEXIPLEX " << VERSION << "\n";
  Rcpp::Rcout << "Setting max flanking sequence edit distance to "
              << max_flank_editdistance << "\n";
  Rcpp::Rcout << "Setting number of threads to " << n_threads << "\n";

  Rcpp::Rcout << "Search pattern:\n";
  for (auto &s : segments) {
    Rcpp::Rcout << s.name << ": " << s.pattern << "\n";
  }
  {
    std::string cb_delim = "";
    Rcpp::Rcout << "CB:Z: tag field: ";
    for (const auto &g : group_map) {
      Rcpp::Rcout << cb_delim << g.first;
      cb_delim = ",";
    }
    for (const auto &s : segments) {
      if (s.type != MATCHED) continue;
      Rcpp::Rcout << cb_delim << s.name;
      cb_delim = ",";
    }
    Rcpp::Rcout << "\n";
  }

  for (const auto &file : {reads_out, stats_out}) {
    if (file_exists(file.get_cstring())) {
      Rcpp::Rcout << "Warning: file " << file.get_cstring()
                  << " already exists, overwriting.\n";
    }
  }

  BGZF *outBgzf = bgzf_open(reads_out.get_cstring(), "w");
  BGZF *statBgzf = bgzf_open(stats_out.get_cstring(), "w");
  bgzf_mt(outBgzf, n_threads, 256);
  bgzf_mt(statBgzf, n_threads, 256);

  // Stats header
  {
    std::string header = "Read";
    for (const auto &s : segments) {
      if (s.type != FIXED)
        header += "\t" + s.name;
    }
    header += "\tFlankEditDist\n";
    write_bgzfstring(statBgzf, header);
  }

  int r_count = 0;
  int r_demultiplexed_count = 0;
  int r_single_match_count = 0;
  int chimeric_count = 0;

  std::unordered_map<std::string, int> barcode_counts;

  int kseq_len = -5;
  for (int i = 0; i < reads_in.size(); i++) {
    Rcpp::Rcout << "Processing file: " << std::string(reads_in(i)) << "\n";

    gzFile gz_reads_in = gzopen(reads_in(i), "r");
    kseq_t *kseq;
    bool is_fastq = false;

    if (!gz_reads_in) {
      Rcpp::stop("Unable to open %s", std::string(reads_in(i)));
    }

    kseq = kseq_init(gz_reads_in);
    kseq_len = kseq_read(kseq);
    if (!(kseq_len >= 0)) {
      Rcpp::stop("Unknown read format");
    }
    is_fastq = (bool)kseq->qual.s;

    gzrewind(gz_reads_in);
    kseq_destroy(kseq);

    kseq = kseq_init(gz_reads_in);

    Rcpp::Rcout << "Searching for barcodes...\n";

    while (kseq_len >= 0) {
      const int buffer_size = 2000;
      std::vector<std::vector<SearchResult>> sr_v(n_threads);
      for (int t = 0; t < n_threads; t++)
        sr_v[t] = std::vector<SearchResult>(buffer_size);

      std::vector<std::thread> threads(n_threads);

      for (int t = 0; t < n_threads; t++) {
        for (int b = 0; b < buffer_size; b++) {
          kseq_len = kseq_read(kseq);
          if (kseq_len < 0) {
            sr_v[t].resize(b);
            for (int t2 = t + 1; t2 < n_threads; t2++)
              sr_v[t2].resize(0);
            if (b > 0) {
              threads[t] = std::thread(search_read, std::ref(sr_v[t]),
                                       std::ref(known_barcodes_map),
                                       std::ref(group_map),
                                       max_flank_editdistance,
                                       std::ref(segments));
            }
            goto print_result;
          }

          SearchResult &sr = sr_v[t][b];
          sr.line = kseq->seq.s;
          sr.read_id = kseq->name.s;
          if (is_fastq)
            sr.qual_scores = kseq->qual.s;

          r_count++;
          if (r_count % 100000 == 0)
            Rcpp::Rcout << r_count / 1000000.0
                        << " million reads processed..\n";
        }

        threads[t] = std::thread(search_read, std::ref(sr_v[t]),
                                 std::ref(known_barcodes_map),
                                 std::ref(group_map), max_flank_editdistance,
                                 std::ref(segments));
      }

    print_result:
      for (int t = 0; t < (int)sr_v.size(); t++) {
        if (!sr_v[t].empty())
          threads[t].join();

        for (int r = 0; r < (int)sr_v[t].size(); r++) {
          auto &res = sr_v[t][r];

          if (res.count > 0)
            r_demultiplexed_count++;
          if (res.count == 1)
            r_single_match_count++;
          if (res.chimeric)
            chimeric_count++;

          // Count first MATCHED segment as the per-read barcode count
          auto count_primary = [&](const Barcode &bc) {
            for (const auto &s : segments) {
              if (s.type == MATCHED) {
                auto it = bc.features.find(s.name);
                if (it != bc.features.end()) {
                  barcode_counts[it->second]++;
                  break;
                }
              }
            }
          };
          for (const auto &bc : res.vec_bc_for)
            count_primary(bc);
          for (const auto &bc : res.vec_bc_rev)
            count_primary(bc);

          print_stats(res.read_id, res.vec_bc_for, statBgzf, segments);
          print_stats(res.read_id, res.vec_bc_rev, statBgzf, segments);

          print_read(res.read_id + "_+", res.line, res.qual_scores,
                     res.vec_bc_for, outBgzf, remove_barcodes, res.chimeric,
                     reverseCompliment, segments, group_map);

          std::reverse(res.qual_scores.begin(), res.qual_scores.end());
          if (remove_barcodes || res.vec_bc_for.empty()) {
            print_read(res.read_id + "_-", res.rev_line, res.qual_scores,
                       res.vec_bc_rev, outBgzf, remove_barcodes, res.chimeric,
                       reverseCompliment, segments, group_map);
          }
        }
      }
    }

    kseq_destroy(kseq);
    gzclose(gz_reads_in);
  }

  bgzf_close(outBgzf);
  bgzf_close(statBgzf);

  Rcpp::Rcout << "Number of reads processed: " << r_count << "\n";
  Rcpp::Rcout << "Number of reads where at least one barcode was found: "
              << r_demultiplexed_count << "\n";
  Rcpp::Rcout << "Number of chimera reads: " << chimeric_count << "\n";
  Rcpp::Rcout << "All done!\n";

  Rcpp::IntegerVector read_counts = Rcpp::IntegerVector::create(
      Rcpp::Named("total reads", r_count),
      Rcpp::Named("demultiplexed reads", r_demultiplexed_count),
      Rcpp::Named("single match reads", r_single_match_count),
      Rcpp::Named("chimera reads", chimeric_count));

  if (kseq_len != -1) {
    const char *msg;
    switch (kseq_len) {
    case -2:
      msg = "truncated quality string";
      break;
    case -3:
      msg = "error reading stream";
      break;
    case -4:
      msg = "overflow error";
      break;
    default:
      msg = "unknown error";
      break;
    }

    Rcpp::stop("Error reading input file %s at read %d: kseq_read returned %d (%s)",
              std::string(reads_in(0)).c_str(), r_count, kseq_len, msg);
  }

  if (barcode_counts.empty())
    return read_counts;

  using pair = std::pair<std::string, int>;
  std::vector<pair> bc_vec;
  bc_vec.reserve(barcode_counts.size());
  for (const auto &kv : barcode_counts)
    bc_vec.push_back(kv);

  std::sort(bc_vec.begin(), bc_vec.end(), [](const pair &l, const pair &r) {
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });

  std::vector<int> hist(bc_vec[0].second);
  std::ofstream out_bc_file;
  out_bc_file.open(bc_out);
  for (auto const &bc_pair : bc_vec) {
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
    if (bc_pair.second > 0)
      hist[bc_pair.second - 1]++;
  }
  out_bc_file.close();

  Rcpp::Rcout << "Reads\tBarcodes\n";
  for (int i = (int)hist.size() - 1; i >= 0; i--)
    if (hist[i] > 0)
      Rcpp::Rcout << i + 1 << "\t" << hist[i] << "\n";

  return read_counts;
}
