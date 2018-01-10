//
// Created by Katie Barr (EI) on 08/01/2018.
//

#include "ComponentPhaser.h"

template<typename T>
void print_vector(std::vector<T> vec){
    for (auto v:vec){
        std::cout << v << " ";
    }
    std::cout << std::endl;
}


ComponentPhaser::ComponentPhaser(SequenceGraph &sg, PairedReadMapper &mapper, std::vector<sgNodeID_t > component, std::vector<std::vector<sgNodeID_t > >  bubbles, MappingParams mapping_params): sg(sg), mapper(mapper) ,component(component), bubbles(bubbles), mapping_params(mapping_params){
    for (int i=0; i < mapper.read_to_tag.size() ; i++){
        if (mapper.read_to_tag[i] == 0){
            std::cout << i << ", "<< " " << mapper.read_to_node[i] << std::endl;;
        }

    }
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
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
    for (int i=0; i < mapper.read_to_tag.size() ; i++){
        if (mapper.read_to_tag[i] == 0){
            std::cout << i << ", " << " " << mapper.read_to_node[i] << std::endl;;
        }

    }
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n 7891 "<< mapper.read_to_tag[7841] << " " << mapper.read_to_tag[8061] << "\n";
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
    int counter = 0;
    for (auto i: mapper.read_to_tag){
        /*if (i == mapper.read_to_tag[0]){
            std::cout << i << ", " << mapper.read_names[counter] << " " << mapper.read_to_node[counter] << std::endl;;
        }
        std::cout << counter << ": " << i << ", ";*/
        counter +=1;

    }
    std::ofstream o1("rtt3.txt");
    o1 << "'";
    for (auto r:mapper.read_to_tag){
        o1 << r << "\n";
    }
    // output file is identical to previous!!! no tags are 0
    o1 << "'" <<std::endl;
    if (mapper.read_to_tag[8061] == 0){std::cout << " yes \n";}
    std::cout << "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n" << mapper.read_to_tag[7841] << " " << mapper.read_to_tag[8061] << " "<< mapper.read_to_tag.size() << " node:  " <<  mapper.read_to_node.size()<< "\n";
    // only care about barcodes for nodes in phaseable bubbles
    for (auto bubble:phaseable_bubbles) {
        print_vector(bubble);
        int mappings_to_bubble = 0;
        int mappings_to_bubble_used = 0;
        std::set<prm10xTag_t> barcodes_bubble_mapped_to;
        int kmer_mappings_to_bubble_used = 0;
        int count = 0;
        for (auto node:bubble) {
            int mappings_to_node_used = 0;

            std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
            std::set<prm10xTag_t> barcodes_node_mapped_to;
            std::cout << mappings.size() << " mappings to node " << node << std::endl;
            mappings_to_bubble += mappings.size();
            // count each barcode mapping to a node in a phaseable bubble
            for (auto mapping:mappings) {
                if (mapping.unique_matches > mapping_params.min_kmer_mappings) {

                        prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
                        if (tag == 0) {std::cout << " id: " << mapping.read_id << " matches " << mapping.unique_matches << " tag:  " << mapper.read_to_tag[mapping.read_id] << " node:  " <<  mapper.read_to_node[mapping.read_id];
                            count +=1;
                            std::ofstream o1("rtt_seroed" + std::to_string(count) + ".txt");
                            o1 << "'";
                            for (auto r:mapper.read_to_tag){
                                o1 << r << "\n";

                            }
                            o1 << "'" <<std::endl;
                        } else {

                            barcode_node_mappings[tag][node] += mapping.unique_matches;
                            mappings_to_node_used += 1;
                            barcodes_mapped_to.insert(tag);
                            barcodes_node_mapped_to.insert(tag);
                            barcodes_bubble_mapped_to.insert(tag);
                            mappings_to_bubble_used += 1;
                            kmer_mappings_to_bubble_used += mapping.unique_matches;
                        }
                    }

            }
            std::cout << std::endl << "node " << mappings_to_node_used << "used " << mappings_to_bubble_used << " mappings, " << kmer_mappings_to_bubble_used << " kmers, " << barcodes_node_mapped_to.size() << " barcodes mapped to \n node in haplotypes:" << std::endl;
            print_vector(node_haplotype_map[node]);
            for (auto h: node_haplotype_map[node]) {
                for (auto tag: barcodes_node_mapped_to) {
                    //std::cout << "tag: " << tag << " h: " << h << " node " << node << std::endl;
                    tags_supporting_haplotypes[h].insert(tag);
                }
            }
        }
        std::cout << barcodes_bubble_mapped_to.size() << " barcodes mapped to bubble" << std::endl;


    }
    for (int i=0; i <possible_haplotypes.size() ; i++) {

        std::cout << i << " " << tags_supporting_haplotypes[i].size() << std::endl;
    }
    // phasing barcodes map to nodes in at least 2 bubbles
    for (auto tag:barcodes_mapped_to){
        std::set<int > bubbles_mapped_to;
        for (auto bubble:phaseable_bubbles) {
            for (auto node:bubble) {
                    bubbles_mapped_to.insert(node_bubble_dict[node]);

            }
        }

        if (bubbles_mapped_to.size() >= 2){
            phasing_barcodes.push_back(tag);
        }
    }
    std::cout << phasing_barcodes.size() << " phasing barcodes of "<< barcodes_mapped_to.size() << std::endl;
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
        if (node_haplotype_map.find(node) != node_haplotype_map.end()){
            for (auto h: node_haplotype_map[node]){
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


void ComponentPhaser::print_haplotype(std::vector<sgNodeID_t > vec){
    for (auto v:vec){
        std::cout << v << " ";
    }
    std::cout << std::endl;

    for (auto v:vec){
        std::cout << sg.oldnames[v] << " ";
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


void ComponentPhaser::print_voting_stats(std::vector<int> vote_totals){

if (!vote_totals.empty()) {
    print_vector(vote_totals);
    int mean = std::accumulate(vote_totals.begin(), vote_totals.end(), 0LL) / vote_totals.size();

    double res = 0;
    for (auto i: vote_totals) {
        res += std::pow(i - mean, 2);
    }
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
    // would be less repetition if used a function for these, but would then iterate over scores 5 times
    // could avoid rhis iterattion by collating these vores above
    // as paired just su, first half to make sure don;t get flippped tie
    for (int i=0; i < scores.size(); i++){
        auto score = scores[i];
        barcode_votes[i] += score.barcode_support;
        overall_votes[i] += score.kmer_support;
        pair_votes[i] += score.pair_support;
        barcode_selecting_votes[i] += score.barcodes_selecting;
    }
    auto barcode_votes_ordered = sort_indexes(barcode_votes);
    auto overall_votes_ordered = sort_indexes(overall_votes);
    auto pair_votes_ordered = sort_indexes(pair_votes);
    auto barcodes_selecting_ordered = sort_indexes(barcode_selecting_votes);
    std::cout << "barcode votex \n";
    print_voting_stats(barcode_votes);
    print_vector(barcode_votes_ordered);
    for (auto ind:barcode_votes_ordered) { std::cout << barcode_votes[ind] << " ";} std::cout << std::endl;
    std::cout << "kmer votex \n";
    print_voting_stats(overall_votes);
    print_vector(overall_votes_ordered);    for (auto ind:overall_votes_ordered) { std::cout << overall_votes[ind] << " ";} std::cout << std::endl;

    std::cout << "barcode pair  votex \n";
    print_voting_stats(pair_votes);
    print_vector(pair_votes_ordered);    for (auto ind:pair_votes_ordered) { std::cout << pair_votes[ind] << " ";} std::cout << std::endl;

    std::cout << "barcodes selecting votes \n";
    print_voting_stats(barcode_selecting_votes);
    print_vector(barcodes_selecting_ordered);    for (auto ind:barcodes_selecting_ordered) { std::cout << barcode_selecting_votes[ind] << " ";} std::cout << std::endl;

    int votes_agreeing = 0;
    // still haven't decded which scores are best....
    if (barcode_votes_ordered[0] == overall_votes_ordered[0]) {
        votes_agreeing += 1;
    }

        if (barcode_votes_ordered[0] == barcodes_selecting_ordered[0]) {
            votes_agreeing += 1;
        }
    if (barcode_votes_ordered[0] == pair_votes_ordered[0]) {
        votes_agreeing += 1;
    }
    if (votes_agreeing >= 2) {
        winning_pair = std::make_pair(barcode_votes_ordered[0],
                                      possible_haplotypes.size() - 1 - barcode_votes_ordered[0]);
        std::cout << votes_agreeing << " scores agree" << std::endl;
        print_haplotype(possible_haplotypes[barcode_votes_ordered[0]]);
        return 1;

    }


};

HaplotypeScore  ComponentPhaser::score_haplotype(size_t index) {
    auto h = possible_haplotypes[index];
    auto pair = possible_haplotypes.size() - index - 1;
    std::cout << "index: "<< index << " pair: " << pair << std::endl;

    auto h_pair = possible_haplotypes[pair];
    HaplotypeScore hs(index);
    hs.barcode_support += tags_supporting_haplotypes[index].size();
    hs.pair_support += tags_supporting_haplotypes[index].size();

    hs.pair_support += tags_supporting_haplotypes[pair].size();

    for (auto tag: tags_supporting_haplotypes[index]) {
        for (auto node: h){
            hs.kmer_support += barcode_node_mappings[tag][node];
            // pair kmer support will always just be all kmers mapping to these contigs!!!

            //hs.pair_kmer_support += barcode_node_mappings[tag][node];

        }
    }

    std::cout << "bs: " << hs.barcode_support << " ps: " << hs.pair_support << " ks: " << hs.kmer_support << std::endl;
     std::cout << std::endl;
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

        }

        std::vector<std::vector<sgNodeID_t >> haps;
        std::cout << P << " combinations to generate from " << N << " bubbles " << std::endl;
        for (size_t m = 0; m < P; m++) {
            // this should hold the index to take from each bubble
            std::vector<size_t> indices(N);
            std::vector<sgNodeID_t> hap;
            size_t m_curr = m;
            for (size_t i = 0; i < N; ++i) {
                indices[i] = m_curr % phaseable_bubbles[i].size();
                hap.push_back(phaseable_bubbles[i][indices[i]]);
                m_curr /= phaseable_bubbles[i].size();
                node_haplotype_map[phaseable_bubbles[i][indices[i]]].push_back(possible_haplotypes.size());
            }

            possible_haplotypes.push_back(hap);
        }
        std::cout << possible_haplotypes.size() << " haplotypes  generated " << std::endl;
        for (auto bubble: phaseable_bubbles) {
            for (auto n: bubble){

                std::cout << "n: " << n << " in haps: ";
                for (auto h: node_haplotype_map[n]) std::cout << h << ", " << mapper.reads_in_node[n].size() << " ";
                std::cout << std::endl<< "next bubble: "<< std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl<< " "<< std::endl;
for (auto h:possible_haplotypes){
    print_haplotype(h);
    for (auto n:h) std::cout  << mapper.reads_in_node[n].size() << " ";
    std::cout << std::endl;
}
        std::cout << "Haplotype nodes size: " << possible_haplotypes.size() << std::endl;
    }


}



    bool ComponentPhaser::node_is_supported(sgNodeID_t node){
    std::vector<ReadMapping> mappings = mapper.reads_in_node[node];
    // determine whether individual node has sufficient mappings to resolve
    // sufficient = mappings with barcodes that map to nodes in other bubbles in this component
    int mappings_with_enough_matches = 0;
        std::map<prm10xTag_t , std::vector< int > > tag_count;
        std::set<prm10xTag_t > tags_used;
    if (mappings.size() > mapping_params.min_node_mappings){
        for (auto mapping:mappings) {
            if (mapping.unique_matches > mapping_params.min_kmer_mappings) {

                // find tag for that mapping
                prm10xTag_t tag = mapper.read_to_tag[mapping.read_id];
                if (mapper.tags_to_nodes[tag].size() > mapping_params.min_nodes_per_tag) {
                    mappings_with_enough_matches += 1;
                    tags_used.insert(tag);

                }
            }
            if (tag_count.find(mapper.read_to_tag[mapping.read_id])!= tag_count.end()) {
                tag_count[mapper.read_to_tag[mapping.read_id]][0] += 1;
                tag_count[mapper.read_to_tag[mapping.read_id]][1] += mapping.unique_matches;
            } else {
                tag_count[mapper.read_to_tag[mapping.read_id]] = {1, mapping.unique_matches};
            }

        }
    }
        if (mappings_with_enough_matches < mapping_params.min_node_mappings_with_enough_matches){
            std::cout << "node: " << node << " mappings: " << mappings.size() << " total tags: " << tag_count.size() << " used: "  << tags_used.size() << std::endl;
            for (auto t: tag_count){
                std::cout << t.first <<": " << t.second[0] << " " << t.second[1] << ",   ";
            }std::cout << std::endl;
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
