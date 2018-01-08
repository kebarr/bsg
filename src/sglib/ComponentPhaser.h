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
    int min_node_mappings_with_enough_matches;
    int min_node_mappings;
    int min_nodes_per_tag;
    int min_kmer_mappings;
};
/*
Barcode::Barcode(static prm10xTag_t barcode) : barcode(barcode){};

class Barcode {
public:
    //Barcode(static prm10xTag_t );
    static prm10xTag_t barcode;
    std::map<sgNodeID_t , int> mappings;
    bool phases_component;
};*/

class ComponentPhaser {
public:
    ComponentPhaser(SequenceGraph &, PairedReadMapper&, std::vector<sgNodeID_t >, std::vector<std::vector<sgNodeID_t > > , MappingParams;

    std::vector<sgNodeID_t > component;
    std::vector<std::vector<sgNodeID_t > > bubbles;

    std::vector<std::vector<sgNodeID_t > > phaseable_bubbles;
    // barcodes which map to nodes from 2 separate bubbles
    std::vector<prm10xTag_t> phasing_barcodes;

    // take component and load nodes that have sufficient barcode support
    void load_bubbles();
    void find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >>, std::map<sgNodeID_t, std::map<prm10xTag_t, int > > , std::map<prm10xTag_t, std::set<sgNodeID_t > > );

private:
    // determine if bubble has enough mappings to be phaseable
    bool bubble_is_supported(std::vector<sgNodeID_t > );
    bool node_is_supported(sgNodeID_t);
    SequenceGraph & sg;
    PairedReadMapper & mapper;
    std::vector<sgNodeID_t > supported_nodes;
    MappingParams mapping_params;
    void load_barcode_mappings();
    std::map<prm10xTag_t, std::map<sgNodeID_t , int>> barcode_node_mappings;
    //std::map<prm10xTag_t, bool> phasing_barcodes;
    std::map<sgNodeID_t , int> node_bubble_dict;

};


#endif //BSG_COMPONENTPHASER_H
