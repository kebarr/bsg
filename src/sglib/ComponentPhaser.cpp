//
// Created by Katie Barr (EI) on 08/01/2018.
//

#include "ComponentPhaser.h"

ComponentPhaser::ComponentPhaser(SequenceGraph &sg, PairedReadMapper &mapper, std::vector<sgNodeID_t > component, std::vector<std::vector<sgNodeID_t > >  bubbles, MappingParams mapping_params): sg(sg), mapper(mapper) ,component(component), bubbles(bubbles), mapping_params(mapping_params){
    for (int i = 0; i < bubbles.size() ; i++){
        for (auto node:bubbles[i]){
            node_bubble_dict[node] = i;
        }
    }
    for (auto node:component){
        if (node_is_supported(node)){
            supported_nodes.push_back(node);
        }
    }
    // order is hard here- node may be supported, but then other bubble mapped to is unsupported
    load_barcode_mappings();
    load_bubbles();
};

// want to reduce number of candidate haplotypes by only using bubbles that mappings can actually phase

void ComponentPhaser::load_barcode_mappings(){
    std::map<prm10xTag_t , int> barcodes_mapped_to;
    // only care about barcodes for nodes in phaseable bubbles
    for (auto bubble:phaseable_bubbles) {
        for (auto node:bubble) {
            std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
            // count each barcode mapping to a node in a phaseable bubble
            for (auto mapping:mappings) {
                if (mapping.unique_matches > mapping_params.min_kmer_mappings) {
                    if (std::find(supported_nodes.begin(), supported_nodes.end(), mapping.node) !=
                        supported_nodes.end()) {
                        prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
                        barcode_node_mappings[tag][node] += mapping.unique_matches;
                    }
                }
            }
        }

    }
    // phasing barcodes map to nodes in at least 2 bubbles
    for (auto tag:barcode_node_mappings){
        std::set<int > bubbles_mapped_to;
        for (auto node:mapper.tags_to_nodes[tag.first]) {
            if (std::find(supported_nodes.begin(), supported_nodes.end(), node) != supported_nodes.end()) {
                bubbles_mapped_to.insert(node_bubble_dict[node]);
            }
        }

        if (bubbles_mapped_to.size() >= 2){
            phasing_barcodes.push_back(tag.first);
        }
    }
};



bool ComponentPhaser::node_is_supported(sgNodeID_t node){
    std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
    // determine whether individual node has sufficient mappings to resolve
    // sufficient = mappings with barcodes that map to nodes in other bubbles in this component
    int mappings_with_enough_matches = 0;
    if (mappings.size() > mapping_params.min_node_mappings){
        for (auto mapping:mappings) {
            if (mapping.unique_matches > mapping_params.min_kmer_mappings) {

                // find tag for that mapping
                prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
                if (mapper.tags_to_nodes[tag].size() > mapping_params.min_nodes_per_tag) {
                    mappings_with_enough_matches += 1;

                }
            }
        }
    }
    return mappings_with_enough_matches > mapping_params.min_node_mappings_with_enough_matches;
};

void ComponentPhaser::load_bubbles(){
    for(auto bubble: bubbles){
        if (bubble_is_supported(bubble)){
            phaseable_bubbles.push_back(bubble);
        }
    }

};

bool ComponentPhaser::bubble_is_supported(std::vector<sgNodeID_t > bubble){
    //bubble is supported if no more thsn 1 node is unsupported

        int support = 0;
        for (auto node: bubble){
            if (std::find(supported_nodes.begin(), supported_nodes.end(), node) != supported_nodes.end()) {
                support +=  1;
            }
        }
        if (support >= bubble.size() -1){

            return true;
        }

};


void ComponentPhaser::find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >>, std::map<sgNodeID_t, std::map<prm10xTag_t, int > > , std::map<prm10xTag_t, std::set<sgNodeID_t > > );
