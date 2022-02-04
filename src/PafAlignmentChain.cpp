#include "PafAlignmentChain.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::to_string;
using std::string;
using std::vector;
using std::sort;
using std::cerr;
using std::cout;
using std::max;

namespace liger2liger {

ChainElement::ChainElement(
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
        bool is_reverse) :
        ref_name(ref_name),
        ref_start(ref_start),
        ref_stop(ref_stop),
        query_start(query_start),
        query_stop(query_stop),
        ref_length(ref_length),
        query_length(query_length),
        residue_matches(residue_matches),
        alignment_length(alignment_length),
        map_quality(mapping_quality),
        is_reverse(is_reverse) {}


uint32_t ChainElement::get_forward_start() const {
    if (is_reverse) {
        return ref_stop;
    } else {
        return ref_start;
    }
}


uint32_t ChainElement::get_forward_stop() const {
    if (is_reverse) {
        return ref_start;
    } else {
        return ref_stop;
    }
}


uint32_t ChainElement::distance_to_end_of_contig() const {
    uint32_t distance;

    if (is_reverse) {
        distance = ref_start;
    } else {
        distance = ref_length - ref_stop;
    }

    return distance;
}


ostream& operator<<(ostream& o, const ChainElement& e) {
    o << '(' << e.query_start << ',' << e.query_stop << ")" << (e.is_reverse ? "-" : "+") << " " << e.ref_name << " "
      << e.ref_start << " " << e.ref_stop << " " << e.ref_length << " " << e.map_quality;
    return o;
}


void AlignmentChains::load_from_paf(path paf_path) {
    ifstream paf_file(paf_path);

    if (not paf_file.good()) {
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string line;

    while (getline(paf_file, line)) {
        add_alignment(line);
    }
}

void print_subchains(
        const PafAlignmentChain& chain,
        const set<pair<size_t, size_t> >& subchain_bounds,
        const string& read_name) {

    cout << "Subchains created for read " << read_name << '\n';

    for (auto& item: subchain_bounds) {
        for (size_t i = item.first; i < item.second; i++) {
            cout << '\t' << chain.chain[i] << '\n';
        }

        cout << '\n';
    }
}


void PafAlignmentChain::add(ChainElement& e) {
    chain.emplace_back(e);
}


bool parse_reversal_token(string& token) {
    bool is_reverse;

    if (token == "-") {
        is_reverse = true;
    } else if (token == "+") {
        is_reverse = false;
    } else {
        throw runtime_error("ERROR: uninterpretable directional symbol is not '-' or '+': " + token);
    }

    return is_reverse;
}

/// Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md
void AlignmentChains::add_alignment(string line) {
    string token;
    string region_name;
    string query_name;
    uint32_t ref_start = 0;
    uint32_t ref_stop = 0;
    uint32_t ref_length = 0;
    uint32_t query_start = 0;
    uint32_t query_stop = 0;
    uint32_t query_length = 0;
    uint32_t map_quality = 0;
    uint32_t residue_matches = 0;
    uint32_t alignment_length = 0;
    uint32_t n_minimizers = 0;
    bool is_reverse = false;

    uint64_t n_delimiters = 0;

    for (char c: line) {
        if (c == '\t') {
            if (n_delimiters == 0) {
                query_name = token;
            } else if (n_delimiters == 1) {
                query_length = stoi(token);
            } else if (n_delimiters == 2) {
                query_start = stoi(token);
            } else if (n_delimiters == 3) {
                query_stop = stoi(token);
            } else if (n_delimiters == 4) {
                is_reverse = parse_reversal_token(token);
            } else if (n_delimiters == 5) {
                region_name = token;
            } else if (n_delimiters == 6) {
                ref_length = stoi(token);
            } else if (n_delimiters == 7) {
                ref_start = stoi(token);
            } else if (n_delimiters == 8) {
                ref_stop = stoi(token);
            } else if (n_delimiters == 9) {
                residue_matches = stoi(token);
            } else if (n_delimiters == 10) {
                alignment_length = stoi(token);
            } else if (n_delimiters == 11) {
                map_quality = stoi(token);
            } else if (n_delimiters == 13) {
                n_minimizers = stoi(token.substr(5, token.size() - 5));

                if (map_quality > min_quality and n_minimizers > min_chain_minimizers) {
                    ChainElement e(
                            region_name,
                            ref_start,
                            ref_stop,
                            query_start,
                            query_stop,
                            ref_length,
                            query_length,
                            residue_matches,
                            alignment_length,
                            map_quality,
                            is_reverse);

                    chains[query_name].add(e);
                }
            }

            token.resize(0);
            n_delimiters++;
        } else if (c == '\n') {
            if (n_delimiters < 13) {
                throw runtime_error("ERROR: file provided does not contain sufficient tab delimiters to be PAF");
            }

            token.resize(0);
            n_delimiters = 0;
        } else {
            token += c;
        }
    }
}


bool compare_chain_elements(ChainElement& a, ChainElement& b) {
    auto midpoint_a = (double(a.query_stop) + double(a.query_start)) / 2;
    auto midpoint_b = (double(b.query_stop) + double(b.query_start)) / 2;

    return midpoint_a < midpoint_b;
}


size_t PafAlignmentChain::size() const {
    return chain.size();
}


void PafAlignmentChain::sort_chain() {
    sort(chain.begin(), chain.end(), compare_chain_elements);
}


uint32_t PafAlignmentChain::compute_distance(ChainElement& a, ChainElement& b) {
    uint32_t distance = 0;
    if (a.ref_name == b.ref_name) {
        auto a_start = a.get_forward_start();
        auto b_start = b.get_forward_start();
        auto a_stop = a.get_forward_stop();
        auto b_stop = b.get_forward_stop();

        // If there is any overlap, set distance to 0
        if ((a_stop > b_start and a_start < b_stop) or (b_stop > a_start and b_start < a_stop)) {
            distance = 0;
        } else {
            distance = abs(int32_t(a_stop) - int32_t(b_start));
        }
    } else {
        // If 2 successive alignments are on different contigs, find the minimum possible distance (+gap penalty)
        int32_t a_to_end = a.distance_to_end_of_contig();
        int32_t b_to_end = b.distance_to_end_of_contig();
        distance = a_to_end + b_to_end + gap_penalty;
    }

    return distance;
}


void PafAlignmentChain::split(set<pair<size_t, size_t> >& subchain_bounds, pair<size_t, size_t> bounds) {
    // For the first recursion, load the result object
    if (subchain_bounds.empty()) {
        bounds = {0, chain.size()};
        subchain_bounds.emplace(bounds);
    }

    size_t start = bounds.first;
    size_t stop = bounds.second;

    uint32_t longest_gap = 0;
    uint32_t gap_index;

    // Iterate and split at largest gap that passes threshold
    // Assume chains have already been sorted by their midpoints
    for (size_t i = start; i < stop - 1; i++) {
        auto gap = compute_distance(chain[i], chain[i + 1]);

        if (gap > longest_gap) {
            longest_gap = gap;
            gap_index = i + 1;
        }
    }

    // Split the bounds if this chain contains a sufficiently large gap
    if (longest_gap > max_gap) {
        subchain_bounds.erase(bounds);

        pair<size_t, size_t> left = {bounds.first, gap_index};
        pair<size_t, size_t> right = {gap_index, bounds.second};

        subchain_bounds.emplace(left);
        subchain_bounds.emplace(right);

        // Recur
        split(subchain_bounds, left);
        split(subchain_bounds, right);
    }
}


void AlignmentChains::split_all_chains() {
    for (auto&[name, chain]: chains) {

        cerr << "before sorting:" << '\n';
        for (auto& item: chain.chain) {
            cerr << item << '\n';
        }
        cerr << '\n';

        // Sort by order of occurrence in query (read) sequence
        chain.sort_chain();

        cerr << "after sorting:" << '\n';
        for (auto& item: chain.chain) {
            cerr << item << '\n';
        }
        cerr << '\n';

        // Do recursive splitting
        set<pair<size_t, size_t> > subchain_bounds;
        chain.split(subchain_bounds);

        print_subchains(chain, subchain_bounds, name);
    }
}

}