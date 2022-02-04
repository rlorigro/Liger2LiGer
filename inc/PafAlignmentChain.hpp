#pragma once

#include "Filesystem.hpp"
#include <ostream>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <list>
#include <map>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::ostream;
using std::vector;
using std::string;
using std::pair;
using std::set;
using std::list;
using std::map;

namespace liger2liger{

class ChainElement {
public:
    string ref_name;
    uint32_t ref_start;
    uint32_t ref_stop;
    uint32_t query_start;
    uint32_t query_stop;
    uint32_t ref_length;
    uint32_t query_length;
    uint32_t residue_matches;
    uint32_t alignment_length;
    uint32_t map_quality;
    bool is_reverse;

    /// Methods ///
    ChainElement(
        string& ref_name,
        uint32_t ref_start,
        uint32_t ref_stop,
        uint32_t query_start,
        uint32_t query_stop,
        uint32_t ref_length,
        uint32_t query_length,
        uint32_t residue_matches,
        uint32_t alignment_length,
        uint32_t mapping_quality,
        bool is_reverse);

    uint32_t distance_to_end_of_contig() const;
    uint32_t get_forward_start() const;
    uint32_t get_forward_stop() const;

    ChainElement()=default;
};


ostream& operator<<(ostream& o, const ChainElement& e);


class PafAlignmentChain {
public:
    vector <ChainElement> chain;

    // Break chains using this threshold for the largest gep between alignments
    static const uint32_t max_gap = 50000;

    // Add this penalty for each time a chain jumps between contigs
    static const uint32_t gap_penalty = 5000;

    /// Methods ///
    PafAlignmentChain()=default;
    void add(ChainElement& e);
    void sort_chain();
    void split(set <pair <size_t, size_t> >& subchain_bounds, pair <size_t, size_t> bounds = {0,0});
    uint32_t compute_distance(ChainElement& a, ChainElement& b);
    size_t size() const;
};


void print_subchains(
        const PafAlignmentChain& chain,
        const set <pair <size_t, size_t> >& subchain_bounds,
        const string& read_name);


class AlignmentChains {
public:
    map <string, PafAlignmentChain> chains;

    // Ignore alignments with mapQ score less than this
    static const uint32_t min_quality = 5;

    // Ignore alignments with a fewer than n minimizers in the chain (using the cm:i:_ tag)
    static const uint32_t min_chain_minimizers = 0;

    /// Methods ///
    AlignmentChains()=default;
    void add_alignment(string line);
    void load_from_paf(path paf_path);
    void split_all_chains();
};

}
