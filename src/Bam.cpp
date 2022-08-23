#include "AlignmentChain.hpp"
#include "Bam.hpp"

#include <stdexcept>
#include <iostream>
#include <vector>

using std::runtime_error;
using std::vector;
using std::cerr;


namespace liger2liger{


Bam::Bam(path bam_path):
    bam_path(bam_path),
    bam_file(nullptr),
//    bam_index(nullptr),
    bam_iterator(nullptr)
{
    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

//    // bam index
//    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == 0) {
//        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
//    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == 0){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    alignment = bam_init1();
}


void Bam::for_alignment_in_bam(const function<void(const string& ref_name, const string& query_name, int32_t query_length, uint8_t map_quality, uint16_t flag)>& f){
    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        string query_name = bam_get_qname(alignment);
        int32_t query_length = alignment->core.l_qseq;

        string ref_name;

        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > 0) {
            ref_name = bam_header->target_name[alignment->core.tid];
        }

        f(ref_name, query_name, query_length, alignment->core.qual, alignment->core.flag);
    }
}


void Bam::for_alignment_in_bam(bool get_cigar, const function<void(SamElement& alignment)>& f){
    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        SamElement e;

        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > -1) {
            e.ref_name = bam_header->target_name[alignment->core.tid];
            e.ref_length = bam_header->target_len[alignment->core.tid];
        }

        e.query_name = bam_get_qname(alignment);
        e.query_length = alignment->core.l_qseq;
        e.mapq = alignment->core.qual;
        e.flag = alignment->core.flag;
        e.ref_start = alignment->core.pos;

        if (get_cigar) {
            auto n_cigar = alignment->core.n_cigar;
            auto cigar_ptr = bam_get_cigar(alignment);
            e.cigars.assign(cigar_ptr, cigar_ptr + n_cigar);
        }

        f(e);
    }
}


//void Bam::for_alignment_in_bam(bool get_cigar, const function<void(ChainElement& alignment)>& f){
//    while (sam_read1(bam_file, bam_header, alignment) >= 0){
//        ChainElement e;
//        e.query_length = alignment->core.l_qseq;
//
//        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
//        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > -1) {
//            e.ref_name = bam_header->target_name[alignment->core.tid];
//        }
//
//        e.ref_start = alignment->core.pos;
//        e.ref_stop = bam_endpos(alignment);
//        e.query_start = bam_q;
//        e.query_stop = ;
//        e.map_quality = alignment->core.qual;
//        e.is_reverse = bam_is_rev(alignment);
//
//        if (get_cigar) {
//            auto n_cigar = alignment->core.n_cigar;
//            auto cigar_ptr = bam_get_cigar(alignment);
//            e.cigars.assign(cigar_ptr, cigar_ptr + n_cigar);
//        }
//
//        f(e);
//    }
//}


Bam::~Bam() {
    hts_close(bam_file);
    bam_hdr_destroy(bam_header);
    bam_destroy1(alignment);
//    hts_idx_destroy(bam_index);
    hts_itr_destroy(bam_iterator);
}


bool Bam::is_first_mate(uint16_t flag){
    return (uint16_t(flag) >> 6) & uint16_t(1);
}


bool Bam::is_second_mate(uint16_t flag){
    return (uint16_t(flag) >> 7) & uint16_t(1);
}


bool Bam::is_not_primary(uint16_t flag){
    return (uint16_t(flag) >> 8) & uint16_t(1);
}


bool Bam::is_primary(uint16_t flag){
    return (not is_not_primary(flag));
}


bool Bam::is_supplementary(uint16_t flag){
    return (uint16_t(flag) >> 11) & uint16_t(1);
}



}
