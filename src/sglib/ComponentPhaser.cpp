//
// Created by Katie Barr (EI) on 08/01/2018.
//

#include "ComponentPhaser.h"

ComponentPhaser::ComponentPhaser(SequenceGraph &sg, PairedReadMapper &mapper, std::vector<sgNodeID_t > component, std::vector<std::vector<sgNodeID_t > >  bubbles, MappingParams mapping_params): sg(sg), mapper(mapper) ,component(component), bubbles(bubbles), mapping_params(mapping_params){

    for (auto bubble:bubbles){
        for (auto node:bubble){
            if (node_is_supported(node)){
                supported_nodes.push_back(node);
            } else {
                unsupported_nodes.push_back(node);
            }
        }
    }
    std::cout << "component of " << component.size() << " nodes  contains " << supported_nodes.size() << " het nodes with supporting mappings \n";

    // order is hard here- node may be supported, but then other bubble mapped to is unsupported
    load_bubbles();
    for (int i = 0; i < phaseable_bubbles.size() ; i++){
        for (auto node:phaseable_bubbles[i]){
            node_bubble_dict[node] = i;
        }
    }
    find_possible_haplotypes();
    std::cout << "mapping data sufficient to phase, in principal, " << phaseable_bubbles.size() << "  bubbles from " << bubbles.size() << " bubbles in total \n";
    load_barcode_mappings();

    for (auto n:component) {
        //hom nodes aren't in bubble, need to check its not an unsupported bubble node
        if (node_bubble_dict.find(n) != node_bubble_dict.end() &
                std::find(unsupported_nodes.begin(), unsupported_nodes.end(), n) != unsupported_nodes.end()) {
            hom_nodes.push_back(n);
        }
    }
};

// want to reduce number of candidate haplotypes by only using bubbles that mappings can actually phase

void ComponentPhaser::load_barcode_mappings(){
    std::set<prm10xTag_t> barcodes_mapped_to;
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
                        for (auto h: bubble_map[node]) {
                            tags_supporting_haplotypes[h].insert(tag);
                        }
                        barcodes_mapped_to.insert(tag);
                    }
                }
            }
        }

    }
    // phasing barcodes map to nodes in at least 2 bubbles
    for (auto tag:barcodes_mapped_to){
        std::set<int > bubbles_mapped_to;
        for (auto node:mapper.tags_to_nodes[tag]) {
            if (std::find(supported_nodes.begin(), supported_nodes.end(), node) != supported_nodes.end()) {
                bubbles_mapped_to.insert(node_bubble_dict[node]);
            }
        }

        if (bubbles_mapped_to.size() >= 2){
            phasing_barcodes.push_back(tag);
        }
    }
};

size_t ComponentPhaser::haplotype_selected_by_barcode(prm10xTag_t barcode){
    // barcode selects a haplotype if it has most mppngs to that haplotype
    std::vector<int> nodes;
    std::vector<int> scores;
    auto nodes_mapped_to = mapper.tags_to_nodes[barcode];
    for (auto node: nodes_mapped_to){
        if (bubble_map.find(node) != bubble_map.end()){
            for (auto h: bubble_map[node]){
                nodes[h] += 1;
                scores[h] += barcode_node_mappings[barcode][node];
            }
        }
    }
    auto highest_scored_haplotype = std::distance(scores, std::max_element(scores.begin(), scores.end());
    auto most_voted_haplotype = std::distance(nodes, std::max_element(nodes.begin(), nodes.end());

};

int ComponentPhaser::phase(){
    auto scores = score_haplotypes();
    
};

std::vector<HaplotypeScore>  ComponentPhaser::score_haplotypes(){
    std::cout << "scoring haplotypes " << std::endl;
    std::vector<HaplotypeScore> scores;
    // kmow number of phasings is aleays even
    for (int index=0; index < possible_haplotypes.size()/2; index++){
    auto h = possible_haplotypes[index];
        auto pair = possible_haplotypes.size() - index - 1;

        auto h_pair = possible_haplotypes[pair];
        //int barcode_support = 0;
        //int kmer_support = 0;
        //int pair_support = 0;
        //int pair_kmer_support = 0;
        HaplotypeScore hs;
        hs.barcode_support += tags_supporting_haplotypes[index].size();
        hs.pair_support += tags_supporting_haplotypes[pair].size();

        for (auto tag: tags_supporting_haplotypes[index]) {
            for (auto node: h){
                hs.kmer_support += barcode_node_mappings[tag][node];

            }
        }
        for (auto tag: tags_supporting_haplotypes[pair]) {

            for (auto node: h_pair) {

                hs.pair_kmer_support += barcode_node_mappings[tag][node];
            }
        }
    scores.push_back(hs);


    }
    for (auto barcode:phasing_barcodes){
        auto winner = haplotype_selected_by_barcode(barcode);
        scores[winner].barcodes_selecting += 1;
    }
    return scores;
};

void ComponentPhaser::find_possible_haplotypes() {
    if (phaseable_bubbles.size() > 1) {
        // algorithm: https://stackoverflow.com/questions/1867030/combinations-of-multiple-vectors-elements-without-repetition
        size_t P = 1;
        auto N = phaseable_bubbles.size();
        for (size_t i = 0; i < N; ++i) {
            P *= phaseable_bubbles[i].size();
            for (auto n: phaseable_bubbles[i]) {
                bubble_map[n].push_back(i);
            }
        }
        std::vector<std::vector<sgNodeID_t >> haps;
        std::cout << P << " combinations to generate from " << N << " bubbles " << std::endl;
        for (size_t m = 0; m < P; m++) {
            // this should hold the index to take from each bubble
            std::vector<size_t> indices(N);
            std::vector<sgNodeID_t> bubble;
            size_t m_curr = m;
            for (size_t i = 0; i < N; ++i) {
                indices[i] = m_curr % phaseable_bubbles[i].size();
                bubble.push_back(phaseable_bubbles[i][indices[i]]);
                m_curr /= phaseable_bubbles[i].size();
            }

            possible_haplotypes.push_back(bubble);
        }
        std::cout << possible_haplotypes.size() << " haplotypes  generated " << std::endl;


        std::cout << "Haplotype nodes size: " << possible_haplotypes.size() << std::endl;
    }


}



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
