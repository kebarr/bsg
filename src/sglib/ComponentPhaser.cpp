//
// Created by Katie Barr (EI) on 08/01/2018.
//

#include "ComponentPhaser.h"

ComponentPhaser::ComponentPhaser(SequenceGraph &sg, PairedReadMapper &mapper, std::vector<sgNodeID_t > component, std::vector<std::vector<sgNodeID_t > >  bubbles, int min_node_mappings=1, int min_nodes_per_tag=2): sg(sg), mapper(mapper) ,component(component), min_node_mappings(min_node_mappings), min_nodes_per_tag(min_nodes_per_tag), bubbles(bubbles){
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

void ComponentPhaser::load_barcode_mappings(){
    std::map<prm10xTag_t , int> barcodes_mapped_to;
    for (auto node:supported_nodes){
        std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
        for (auto mapping:mappings){
            prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
            barcode_node_mappings[tag][node] += 1;
        }
    }

    for (auto tag:barcode_node_mappings){
        int nodes_mapped_to = 0;
        for (auto node:mapper.tags_to_nodes[tag.first]) {
            if (std::find(supported_nodes.begin(), supported_nodes.end(), node) != supported_nodes.end()) {
                nodes_mapped_to += 1;
            }
        }
        // shouldn't need this as shoul
        if (nodes_mapped_to > min_nodes_per_tag){
            phasing_barcodes.push_back(tag.first);
        }
    }
};



bool ComponentPhaser::node_is_supported(sgNodeID_t node){
    std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
    // determine whether individual node has sufficient mappings to resolve
    // sufficient = mappings with barcodes that map to nodes in other bubbles in this component
    if (mappings.size() > min_node_mappings){
        for (auto mapping:mappings) {
            // find tag for that mapping
            prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
            if (mapper.tags_to_nodes[tag].size() > min_nodes_per_tag) {
                for (auto other_node: mapper.tags_to_nodes[tag]) {
                        // tag must be mapped to node in other  bubble of component to be useful
                        if (std::find(component.begin(), component.end(), node) != component.end()) {
                            int bubble = node_bubble_dict[node];

                            if (bubble != node_bubble_dict[other_node]) {

                                return true;
                        }
                    }
                }

            }
        }
    }
    return  false;
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
