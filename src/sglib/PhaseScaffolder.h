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
    // max bubbles per phase block - too many explodes cobinations, and gets very weak signal as barcodes don't map across entire component
    // min_barcodes_mapping is how many votes are required to decide a phasing
    void phase_components(int max_bubbles=12, int min_barcodes_mapping=10);
    void load_mappings(std::string , std::string , std::string, uint64_t , std::string);
    void load_mappings_from_file(std::string );

        PairedReadMapper mapper;
    void intersect_phasings();


        private:
    void print_barcode_stats();
    std::map<prm10xTag_t, std::set<sgNodeID_t > > barcode_node_mappings_int;
    int phase_component (std::vector<std::vector<sgNodeID_t >>, std::vector<sgNodeID_t >, int min_barcodes_mapping=2);
    std::map< sgNodeID_t, std::map< prm10xTag_t, int > > node_tag_mappings;
    void sum_node_tag_mappings(std::vector< std::vector<prm10xTag_t> >, int min_tag_count=1);
    std::vector<HaplotypeScorer> phased_components;
    std::vector<HaplotypeScorer> partial_phased_components;
    std::vector<std::vector<std::vector<sgNodeID_t > > > split_component(std::vector<std::vector<sgNodeID_t >>, int max_bubbles=12);


    };


#endif //SG_PHASESCAFFOLDER_H
