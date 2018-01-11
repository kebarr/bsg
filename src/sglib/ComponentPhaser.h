//
// Created by Katie Barr (EI) on 08/01/2018.
//

#ifndef BSG_COMPONENTPHASER_H
#define BSG_COMPONENTPHASER_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <numeric>

#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>

struct MappingParams{
    // arbitrary defaults
    int min_node_mappings_with_enough_matches = 2;
    int min_node_mappings = 2;
    int min_nodes_per_tag = 2;
    int min_kmer_mappings = 2;
};

struct HaplotypeScore{
    HaplotypeScore(size_t index) : index(index){};
    HaplotypeScore(){};
    int barcode_support = 0;
    int kmer_support = 0;
    int barcodes_against_pair= 0;
    int pair_support = 0;
    //int pair_kmer_support = 0;
    int barcodes_selecting = 0;
    size_t  index;

};


class ComponentPhaser {
public:
    ComponentPhaser(SequenceGraph &, PairedReadMapper&, std::vector<sgNodeID_t >, std::vector<std::vector<sgNodeID_t > > , MappingParams);

    std::vector<sgNodeID_t > component;
    std::vector<std::vector<sgNodeID_t > > bubbles;
    std::vector<sgNodeID_t > hom_nodes;

    std::vector<std::vector<sgNodeID_t > > phaseable_bubbles;
    std::vector<std::vector<sgNodeID_t > > possible_haplotypes;
    // barcodes which map to nodes from 2 separate bubbles
    std::vector<prm10xTag_t> phasing_barcodes;

    // take component and load nodes that have sufficient barcode support
    void load_bubbles();
    size_t phase();
    void find_possible_haplotypes();
    std::map<size_t , std::set<prm10xTag_t >> tags_supporting_haplotypes;
    std::pair<size_t , size_t > winning_pair;

private:
    std::vector<HaplotypeScore>  score_haplotypes();
    std::map<sgNodeID_t , std::vector<size_t > > node_haplotype_map;
    // determine if bubble has enough mappings to be phaseable
    bool bubble_is_supported(std::vector<sgNodeID_t > );
    bool node_is_supported(sgNodeID_t);
    SequenceGraph & sg;
    PairedReadMapper & mapper;
    std::vector<sgNodeID_t > supported_nodes;
    std::vector<sgNodeID_t > unsupported_nodes;
    HaplotypeScore  score_haplotype(size_t);

    void print_voting_stats(std::vector<int>);

                                             MappingParams mapping_params;
    void load_barcode_mappings();
     std::map<prm10xTag_t,std::map<sgNodeID_t , int>> barcode_node_mappings;
    //std::map<prm10xTag_t, bool> phasing_barcodes;
    std::map<sgNodeID_t , int> node_bubble_dict;
    std::vector<HaplotypeScore> haplotype_scores;

    std::vector<size_t> haplotype_selected_by_barcode(prm10xTag_t );

    void print_haplotype(std::vector<sgNodeID_t > vec);
};


#endif //BSG_COMPONENTPHASER_H
