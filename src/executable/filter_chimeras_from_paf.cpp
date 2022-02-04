#include "PafAlignmentChain.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <queue>
#include <cmath>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::queue;
using std::pair;
using std::cerr;
using std::cout;
using std::abs;

using liger2liger::AlignmentChains;
using liger2liger::PafAlignmentChain;
using liger2liger::ChainElement;


bool chain_is_palindromic(const PafAlignmentChain& chain, const pair <size_t, size_t>& bounds){
    bool is_palindromic = false;

    uint32_t prev_midpoint = 0;
    bool prev_reversal;
    uint32_t n_strand_reversals = 0;
    size_t reversal_index = 0;

    uint32_t left_length = 0;
    uint32_t right_length = 0;

    set<string> left_contigs;
    set<string> right_contigs;

    // Iterate and accumulate chains for continuous strands. Chain must be sorted in order of query coordinates.
    for (size_t i=bounds.first; i < bounds.second; i++) {
        const ChainElement& c = chain.chain[i];

        if (i == 0){
            prev_reversal = c.is_reverse;
        }

        auto midpoint = (double(c.query_stop) - double(c.query_stop))/2;

        if (midpoint < prev_midpoint){
            throw runtime_error("ERROR: attempting to classify palindromic chain for unsorted chain");
        }

        prev_midpoint = midpoint;

        // New strands start whenever the "reversal" flag flips
        if (prev_reversal != c.is_reverse){
            n_strand_reversals++;

            if (n_strand_reversals == 0){
                reversal_index = i;
            }
        }

        if (n_strand_reversals == 0){
//            left_length += c.query_stop - c.query_start;
            left_length += c.alignment_length;
            left_contigs.emplace(c.ref_name);
        }
        if (n_strand_reversals == 1){
//            right_length += c.query_stop - c.query_start;
            right_length += c.alignment_length;
            right_contigs.emplace(c.ref_name);
        }

    }

    // A palindrome should have only 1 strand reversal and both strands should map to one contig
    if (n_strand_reversals == 1 and left_contigs.size() == 1 and right_contigs.size() == 1){

        // Both strands should map to THE SAME contig
        if (*left_contigs.begin() == *right_contigs.begin()){
            double length_ratio = double(left_length)/double(right_length);

            cerr << '\t' << "length ratio: " << length_ratio << '\n';

        }
    }

    return is_palindromic;
}


void filter_paf(path paf_path){
    AlignmentChains alignment_chains;
    alignment_chains.load_from_paf(paf_path);

    path chimer_id_path = paf_path;
    chimer_id_path.replace_extension("chimeric_reads.txt");
    ofstream chimer_id_file(chimer_id_path);

    path non_chimer_id_path = paf_path;
    non_chimer_id_path.replace_extension("non_chimeric_reads.txt");
    ofstream non_chimer_id_file(non_chimer_id_path);

    cerr << "Writing chimeric reads to file: " << chimer_id_path << '\n';
    cerr << "Writing non-chimeric reads to file: " << non_chimer_id_path << '\n';

    vector<size_t> non_chimer_lengths;
    vector<size_t> chimer_lengths;

    path non_chimer_lengths_path = paf_path;
    path chimer_lengths_path = paf_path;
    path chimer_subchains_lengths_path = paf_path;
    path chimer_subchains_path = paf_path;

    non_chimer_lengths_path.replace_extension("non_chimer_lengths.txt");
    chimer_lengths_path.replace_extension("chimer_lengths.txt");
    chimer_subchains_lengths_path.replace_extension("chimer_subchains_lengths.txt");
    chimer_subchains_path.replace_extension("chimer_subchains.txt");

    ofstream non_chimer_lengths_file(non_chimer_lengths_path);
    ofstream chimer_lengths_file(chimer_lengths_path);
    ofstream chimer_subchains_lengths_file(chimer_subchains_lengths_path);
    ofstream chimer_subchains_file(chimer_subchains_path);

    for (auto& [name, chain]: alignment_chains.chains) {
        // Sort by order of occurrence in query (read) sequence
        chain.sort_chain();

        // Do recursive splitting to find the index bounds of sub-chains
        set <pair <size_t, size_t> > subchain_bounds;
        chain.split(subchain_bounds);

        if (subchain_bounds.size() > 1) {
            vector<uint32_t> chain_lengths;
            chain_lengths.resize(subchain_bounds.size());

            // Iterate subchains created by splitting
            for (auto &item: subchain_bounds) {
                for (uint32_t i = item.first; i < item.second; i++) {
                    uint32_t length = abs(int32_t(chain.chain[i].query_stop) - int32_t(chain.chain[i].query_start));
                    chimer_subchains_lengths_file << length << '\n';
                }
            }

            chimer_lengths_file << chain.chain[0].query_length << '\n';

            chimer_id_file << name << '\n';

            chimer_subchains_file << name << '\t';
            for (auto& item: subchain_bounds) {
                chimer_subchains_file << '(' << chain.chain[item.first].query_start << ',' << chain.chain[item.second - 1].query_stop << ")" << ',';
            }
            chimer_subchains_file << '\n';

//            print_subchains(chain, subchain_bounds, name);
        }
        else{
            non_chimer_id_file << name << '\n';
            non_chimer_lengths_file << chain.chain[0].query_length << '\n';
        }
    }
}


int main(int argc, char* argv[]){
    path paf_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--paf_path",
            paf_path,
            "File path of PAF file containing alignments to some reference")
            ->required();


    filter_paf(paf_path);

    return 0;
}
