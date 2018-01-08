//
// Created by Katie Barr (EI) on 08/01/2018.
//

#include "ComponentPhaser.h"

ComponentPhaser::ComponentPhaser(SequenceGraph &sg, PairedReadMapper &mapper, std::vector<sgNodeID_t > component, int min_node_mappings=1, int min_nodes_per_tag=2): sg(sg), mapper(mapper) ,component(component), min_node_mappings(min_node_mappings), min_nodes_per_tag(min_tags_per_node){
    for (auto node:component){
        if (node_is_supported(node)){
            supported_nodes.push_back(node);
        }
    }
};

void ComponentPhaser::load_barcode_mappings(){
    std::set<prm10xTag_t > barcodes_mapped_to;
    for (auto node:supported_nodes){
        std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
        for (auto mapping:mappings){
            prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
            barcodes_mapped_to.insert(tag);
        }
    }

    for (auto tag:barcodes_mapped_to){
        int nodes_mapped_to = 0;
        for (auto node:mapper.tags_to_nodes[tag]) {
            if (std::find(supported_nodes.begin(), supported_nodes.end(), node) != supported_nodes.end()) {
                barcode_node_mappings[tag][node] += 1;
                nodes_mapped_to += 1;
            }
        }
        // shouldn't need this as shoul
        if (nodes_mapped_to > min_nodes_per_tag){
            phasing_barcodes.push_back(tag);
        }
    }
};



bool ComponentPhaser::node_is_supported(sgNodeID_t node){
    std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
    // determine whether individual node has sufficient mappings to resolve
    //
    if (mappings.size() > min_node_mappings){
        for (auto mapping:mappings) {
            // find tag for that mapping
            prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
            if (mapper.tags_to_nodes[tag].size() > min_nodes_per_tag) {
                for (auto node: mapper.tags_to_nodes[tag]) {
                    // tag must be mapped to other node in component to be useful
                    if (std::find(component.begin(), component.end(), node) != component.end()){
                        return true;
                    }
                }

            }
        }
    }
    return  false;
};

void ComponentPhaser::load_bubbles(std::vector<std::vector<sgNodeID_t >> bubbles){
    for(auto bubble: bubbles){
        if (bubble_is_supported(bubble)){
            phaseable_bubbles.push_back(bubble);
        }
    }
};


void ComponentPhaser::find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >>, std::map<sgNodeID_t, std::map<prm10xTag_t, int > > , std::map<prm10xTag_t, std::set<sgNodeID_t > > );
