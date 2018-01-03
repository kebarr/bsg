//
// Created by Katie Barr (EI) on 23/11/2017.
//

#ifndef SG_PHASESCAFFOLDER_H
#define SG_PHASESCAFFOLDER_H

#include <string>


#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>
#include <sglib/HaplotypeScorer.hpp>
/**
 *
 *
 * @brief
 * Find haplotypes for each component and join into phased blocks
 *
 */

class PhaseScaffolder {
public:
    PhaseScaffolder(SequenceGraph &);
    SequenceGraph & sg;

    void output_bubbles(std::string);
    void phase_components(int max_bubbles=12);
    void load_mappings(std::string , std::string , std::string, uint64_t , std::string);
    void load_mappings_from_file(std::string );

        PairedReadMapper mapper;
    void intersect_phasings();


        private:
    void print_barcode_stats();
    int phase_component (std::vector<std::vector<sgNodeID_t >>);
    std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings;
    void sum_node_tag_mappings(std::vector< std::vector<prm10xTag_t> >, int min_tag_count=1);
    std::vector<HaplotypeScorer> phased_components;
    std::vector<HaplotypeScorer> partial_phased_components;

};


#endif //SG_PHASESCAFFOLDER_H
