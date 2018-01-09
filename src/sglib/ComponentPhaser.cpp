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

//from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

    return idx;
}

size_t ComponentPhaser::haplotype_selected_by_barcode(prm10xTag_t barcode){
    // barcode selects a haplotype if it has most mppngs to that haplotype
    std::vector<int> nodes(possible_haplotypes.size());
    std::vector<int> scores(possible_haplotypes.size());
    auto nodes_mapped_to = mapper.tags_to_nodes[barcode];
    for (auto node: nodes_mapped_to){
        if (bubble_map.find(node) != bubble_map.end()){
            for (auto h: bubble_map[node]){
                nodes[h] += 1;
                scores[h] += barcode_node_mappings[barcode][node];
            }
        }
    }
    auto scores_ordered = sort_indexes(scores);
    auto nodes_ordered = sort_indexes(nodes);
    auto highest_scored_haplotype = scores_ordered[0];
    auto most_voted_haplotype = nodes_ordered[0];
    if (highest_scored_haplotype == most_voted_haplotype){
        return highest_scored_haplotype;
    } else if (highest_scored_haplotype = possible_haplotypes.size() - 1 - most_voted_haplotype){
        return  highest_scored_haplotype;
    } else {
        return  possible_haplotypes.size();
    }

};

void print_vector(std::vector<int> vec){
    for (auto v:vec){
        std::cout << v << " ";
    }
    std::cout << std::endl;
}


void print_haplotype(std::vector<sgNodeID_t > vec){
    for (auto v:vec){
        std::cout << v << " ";
    }
    std::cout << std::endl;
}


double avg(std::vector<int> v){
    if (v.size() > 0) {
        return std::accumulate(v.begin(), v.end(), 0LL) / v.size();
    }
    return 0.0;
}

double stdev(std::vector<int> v, double mean){
    if (v.size() > 0) {
        double res = 0;
        for (auto i: v) {
            res += std::pow(i - mean, 2);
        }
        return std::pow(res / v.size(), 0.5);
    }
    return 0.0;
}
float calculateSD(std::vector<int> data)
{
    int sum = 0;
    float mean, standardDeviation = 0.0;

    int i;

    for(i = 0; i < data.size(); ++i)
    {
        sum += data[i];
        std::cout << sum << " i: " << data[i]<<", ";
    }
    std::cout << std::endl;

    mean = sum/data.size();

    for(i = 0; i < data.size(); ++i){
        standardDeviation += pow(data[i] - mean, 2);
    std::cout << standardDeviation << " i: " << data[i]<<", " << pow(data[i] - mean, 2) <<", " ;}
    std::cout << std::endl << standardDeviation  << std::endl ;

    return sqrt(standardDeviation / data.size());
}


void ComponentPhaser::print_voting_stats(std::vector<int> vote_totals){

if (!vote_totals.empty()) {
    print_vector(vote_totals);
    int mean = std::accumulate(vote_totals.begin(), vote_totals.end(), 0LL) / vote_totals.size();

    std::cout << " mean" << mean << std::endl;
    std::cout << calculateSD(vote_totals) << std::endl;;
    double res = 0;
    for (auto i: vote_totals) {
        std::cout << res << " i: " << i <<", " << std::pow(i - mean, 2)  <<", "  << i - mean <<", " << mean <<", ";
        res += std::pow(i - mean, 2);
    }
    std::cout << std::endl << res << std::endl ;
    auto stdev = std::pow(res / vote_totals.size(), 0.5);
    auto support_max = std::max_element(vote_totals.begin(), vote_totals.end());
    auto support_min = std::min_element(vote_totals.begin(), vote_totals.end());
    std::cout << "Max: " << *support_max << " min: " << *support_min << " mean: " << mean << " stdev: " << stdev <<" for " << vote_totals.size() << " haplotypes " <<std::endl;
    }
}


size_t ComponentPhaser::phase(){
    auto scores = score_haplotypes();
    std::vector<int> barcode_votes(possible_haplotypes.size());
    std::vector<int> overall_votes(possible_haplotypes.size());
    std::vector<int> pair_votes(possible_haplotypes.size());

    std::vector<int> barcode_selecting_votes(possible_haplotypes.size());
    std::vector<int> pair_overall_votes(possible_haplotypes.size());
    // would be less repetition if used a function for these, but would then iterate over scores 5 times
    // could avoid rhis iterattion by collating these vores above
    for (int i=0; i < scores.size(); i++){
        auto score = scores[i];
        barcode_votes[i] += score.barcodes_selecting;
        overall_votes[i] += score.kmer_support;
        pair_votes[i] += score.pair_support;
        pair_overall_votes[i] += score.pair_kmer_support;
        barcode_selecting_votes[i] += score.barcodes_selecting;
    }
    auto barcode_votes_ordered = sort_indexes(barcode_votes);
    auto overall_votes_ordered = sort_indexes(overall_votes);
    auto pair_votes_ordered = sort_indexes(pair_votes);
    auto pair_overall_votes_ordered = sort_indexes(pair_overall_votes);
    auto barcodes_selecting_ordered = sort_indexes(barcode_selecting_votes);
    std::cout << "barcode votex \n";
    print_voting_stats(barcode_votes);
    std::cout << "kmer votex \n";
    print_voting_stats(overall_votes);
    std::cout << "barcode pair  votex \n";
    print_voting_stats(pair_votes);
    std::cout << "pair kmer votex \n";
    print_voting_stats(pair_overall_votes);
    int votes_agreeing = 0;
    // still haven't decded which scores are best....
    if (barcode_votes_ordered[0] == overall_votes_ordered[0]) {
        votes_agreeing += 1;
    }

        if (barcode_votes_ordered[0] == barcodes_selecting_ordered[0]) {
            votes_agreeing += 1;
        }
    if (barcode_votes_ordered[0] == pair_overall_votes_ordered[0]) {
        votes_agreeing += 1;
    }
    if (barcode_votes_ordered[0] == pair_votes_ordered[0]) {
        votes_agreeing += 1;
    }
    if (votes_agreeing >= 3) {
        winning_pair = std::make_pair(barcode_votes_ordered[0],
                                      possible_haplotypes.size() - 1 - barcode_votes_ordered[0]);
        std::cout << votes_agreeing << " scores agree" << std::endl;
        print_haplotype(possible_haplotypes[barcode_votes_ordered[0]]);
        return barcode_votes_ordered[0];

    }


};

HaplotypeScore  ComponentPhaser::score_haplotype(size_t index) {
    auto h = possible_haplotypes[index];
    auto pair = possible_haplotypes.size() - index - 1;

    auto h_pair = possible_haplotypes[pair];
    HaplotypeScore hs(index);
    hs.barcode_support += tags_supporting_haplotypes[index].size();
    hs.pair_support += tags_supporting_haplotypes[index].size();

    hs.pair_support += tags_supporting_haplotypes[pair].size();

    for (auto tag: tags_supporting_haplotypes[index]) {
        for (auto node: h){
            hs.kmer_support += barcode_node_mappings[tag][node];

        }
    }
    for (auto tag: tags_supporting_haplotypes[pair]) {

        for (auto node: h_pair) {

            hs.pair_kmer_support += barcode_node_mappings[tag][node];
            hs.kmer_support += barcode_node_mappings[tag][node];
        }
    }
    return hs;
}

std::vector<HaplotypeScore>  ComponentPhaser::score_haplotypes(){
    std::cout << "scoring haplotypes " << std::endl;
    std::vector<HaplotypeScore> scores(possible_haplotypes.size());
    // kmow number of phasings is aleays even
    for (int index=0; index < possible_haplotypes.size()/2; index++){
        auto pair = possible_haplotypes.size() - index - 1;
    auto hs = score_haplotype(index);
       auto hs_pair=  score_haplotype(pair);
    scores[index] = hs;
        scores[pair] = hs_pair;


    }
    for (auto barcode:phasing_barcodes){
        auto winner = haplotype_selected_by_barcode(barcode);
        if (winner != possible_haplotypes.size()) {
            scores[winner].barcodes_selecting += 1;
        }
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
