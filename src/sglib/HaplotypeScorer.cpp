//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "HaplotypeScorer.hpp"

void print_vector(std::vector<std::string> vec){
    int count = 0;
    for (auto a: vec){
        std::cout << count << ": " << a << " ";
        count++;
    }
    std::cout << std::endl;
}

void print_id_vector(std::vector<sgNodeID_t > vec){
    int count = 0;
    for (auto a: vec){
        std::cout << count << ": " << a << " ";
        count++;
    }
    std::cout << std::endl;
}

void print_int_vector(std::vector<int> vec){
    int count = 0;
    for (auto a: vec){
        std::cout << count << ": " << a << " ";
        count++;
    }
    std::cout << std::endl;
}

void print_pair_int_vector(std::vector<std::pair<int, int> > vec){
    for (auto a: vec){
        std::cout << std::get<0>(a) << " " << std::get<1>(a);
    }
    std::cout << std::endl;
}


void print_pair_int_map(std::map<std::pair<int, int> , int> map){
    for (auto a: map){
        std::cout << std::get<0>(a.first) << " " << std::get<1>(a.first) << ": " << a.second << " ";
    }
    std::cout << std::endl;
}

void HaplotypeScorer::remove_nodes_with_no_barcode_support(std::vector< std::vector<prm10xTag_t> > tags_in_node){
    std::cout<< " updating haplotypes to include only supported nodes"<<std::endl;

    // ONLY COUNT NODES THAT ARE VOTED FOR- EG H1 IS 1,2,3,4 BUT NO VOTES FOR 4, AND H2 IS 1,2,3,5 WITH NONE FOR 5- THESE ARE SAME!
    std::set<std::set<sgNodeID_t > > haplotype_ids_remove_unmapped_nodes;
    for (int i=0; i< haplotype_ids.size() ; i++){
        auto h= haplotype_ids[i];
        std::set<sgNodeID_t >  supported_nodes;
        for (auto n:h){
            if (!tags_in_node[n].empty()){
                supported_nodes.insert(n);
            }
        }
        //std::cout << "supported nodes : " << supported_nodes.size() <<std::endl;
        haplotype_ids_remove_unmapped_nodes.insert(supported_nodes);
    }
    // removing unsupported nodes should collapse combinations which don't have read mappinhgs for each node
    // this shouls stop haplotypes from drawing
    haplotype_ids.clear();
    for (auto h: haplotype_ids_remove_unmapped_nodes){
        std::vector<sgNodeID_t > haplotype;
        for (auto n:h){
            haplotype.push_back(n);
        }
        haplotype_ids.push_back(haplotype);
    }
    std::cout<< "now " << haplotype_ids.size() << " haplotypes"<<std::endl;
    this-> haplotype_ids = haplotype_ids;
}

/**
 *
 *
 * @brief
 * sum each barcode's score for each haplotype
 * @return
 * Returns number of barcodes that can be used to phase this component- i.e. map to more than 2 het sites
 * For each barcode, calculate how many nodes it shares with each haplotype- i.e. how many nodes it ,aps to contained in each haplotype
 * then for each node in each haplotype, find barcodes mapping
 */

int HaplotypeScorer::decide_barcode_haplotype_support(std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings, std::map<prm10xTag_t, std::vector<sgNodeID_t > > barcode_node_mappings){

    int total_mappings = 0;
    int total_barcodes = 0;
    std::map<size_t , std::map< prm10xTag_t, int>> barcode_haplotype_shared;// number of nodes that a barcode maps to in each hap
    size_t shared_nodes = 0;

    // for each barcode, calculate which haps it shares nodes with
    for (auto m: barcode_node_mappings){
        auto barcode = m.first;
        auto nodes = m.second;
        for (int i = 0; i < haplotype_ids.size() ; i++) {
            for (auto l: nodes) {
                if (std::find(haplotype_ids[i].begin(), haplotype_ids[i].end(), l) != haplotype_ids[i].end()) {
                    barcode_haplotype_shared[i][barcode] += 1;

                }
            }

        }
        }

    haplotype_barcodes_total_mappings.resize(haplotype_ids.size());
    haplotype_barcodes_supporting.resize(haplotype_ids.size());
    // previously found all nodes that a barcode maps to and only considered ones which mapped to enough of the haploype
    // for each haplotype, loop over each node and sum support
    for (int i = 0; i < haplotype_ids.size() ; i++){
        // to be same as previous, need to check each barcode counted maps to enough nodes in haplotype
        for (auto node: haplotype_ids[i]){
            if (node_tag_mappings[node].size() > 1) {// if more than 1 barcode maps tp this node
                //node_tag_mappings is a map: node_id:barcode_id:number of mappings from that barcode to that node
                for (auto f:node_tag_mappings[node]) {// assign each barcode mapped to this node to the relevant haplotype
                    // if this barcode maps to enough nodes in this haplotype
                    // f.first is barcode
                    auto barcode = f.first;
                    if (barcode_haplotype_shared[i][barcode] > node_tag_mappings[node][barcode]/2) {
                        haplotype_barcodes_supporting[i] += 1;

                            haplotype_barcodes_total_mappings[i] += f.second;
                            total_mappings += f.second;
                        total_barcodes += 1;
                        if (barcodes_supporting_haplotype[i].size() > 0){
                            barcodes_supporting_haplotype[i].push_back(barcode);
                            } else {
                            barcodes_supporting_haplotype[i]= {barcode};

                        }

                    }
                }
            }
        }
    }
for (int i=0; i < haplotype_barcodes_supporting.size() ; i++){
    if (haplotype_barcodes_supporting[i] > 0 || haplotype_barcodes_total_mappings[i] > 0) {
        std::cout << "i: " << i << " votes " << haplotype_barcodes_supporting[i] << " kmers "
                  << haplotype_barcodes_total_mappings[i] << std::endl;
        for (auto h: haplotype_ids[i]) {
            std::cout << h << " ";
        }
        std::cout << std::endl;
    }
}
    this->haplotype_barcodes_supporting = haplotype_barcodes_supporting;
    this->haplotype_barcodes_total_mappings = haplotype_barcodes_total_mappings;
    this->barcodes_supporting_haplotype = barcodes_supporting_haplotype;
    std::cout << haplotype_barcodes_supporting[0] << " " << haplotype_barcodes_supporting[1] << std::endl;
    std::cout << haplotype_barcodes_total_mappings[0] << " " << haplotype_barcodes_total_mappings[1] << std::endl;

              std::cout << "haplotype_barcodes_supporting.size() "<< haplotype_barcodes_supporting.size() <<  std::endl;
    std::cout << "haplotype_barcodes_total_mappings.size() "<< haplotype_barcodes_total_mappings.size() <<  std::endl;

    std::cout << "barcodes_supporting_haplotype.size() "<< barcodes_supporting_haplotype.size() <<  std::endl;

    std::cout << "Calculated haplotype support for each barcode, total mappings: " << total_mappings<< "total barcodes: " << total_barcodes <<  std::endl;
    return  total_barcodes;
}


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

void HaplotypeScorer::print_voting_stats(std::vector<int> vote_totals){
    if (vote_totals.size() > 0) {
        auto mean = std::accumulate(vote_totals.begin(), vote_totals.end(), 0LL) / vote_totals.size();

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
 /*
 *
 * @brief
 * Decide winning phasing
 * @return
 * Calculate, for each haplotype, a number of metrics that represent how well the 10x reads support that phasing
 * loop over each barcode in turn, then each phasing pair
 * \todo decide criteria for winning haploype
 *
 */
 int HaplotypeScorer::score_haplotypes(std::vector<std::string> oldnames){
     std::cout << "scoring " << haplotype_ids.size() << " possible phasings \n";
     // need pairs with max barcodes + kmers support
     auto ordered_haplotype_barcodes_supporting = sort_indexes(haplotype_barcodes_supporting);

     auto ordered_haplotype_barcodes_total_mappings = sort_indexes(haplotype_barcodes_total_mappings);
     std::cout << ordered_haplotype_barcodes_supporting[0] << " won with score " << haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[0]] << " " << ordered_haplotype_barcodes_supporting[1] << " runner up with score " << haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[1]] << std::endl;
     std::cout << ordered_haplotype_barcodes_total_mappings[0] << " won with score "  << haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[0]] << " " << ordered_haplotype_barcodes_total_mappings[1] << " runner up  with score "  <<haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[1]] << std::endl;

     /*for (int i = 0 ; i < ordered_haplotype_barcodes_supporting.size() ; i++){

         if (haplotype_barcodes_supporting[i] > 0){
             std::cout << "h: " << haplotype_barcodes_supporting[i]  << " " << i << " ";
         }

     }
     std::cout << std::endl;
     for (int i = 0 ; i < ordered_haplotype_barcodes_total_mappings.size() ; i++){

         if (haplotype_barcodes_total_mappings[i] > 0){
             std::cout << "h: " << haplotype_barcodes_total_mappings[i]  << " " << i << " ";
         }
     }
     std::cout << std::endl;*/
    // if there was a vote (rather than 0 barcodes for each haplotype) and no tie, check winner
     if (haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[0]] >0 && haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[0]] > 0
         //&& haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[0]] != haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[1]] &&
           //  haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[0]] != haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[1]]
             ){
         std::cout << "barcode support stats: " << std::endl;
         print_voting_stats(haplotype_barcodes_supporting);

         std::cout << "kmer support stats: " << std::endl;
         print_voting_stats(haplotype_barcodes_total_mappings);
         auto winner = ordered_haplotype_barcodes_supporting[0];
         auto winner_pair = haplotype_ids.size() - 1 -
                            ordered_haplotype_barcodes_supporting[0];
         auto winner_mappings = ordered_haplotype_barcodes_total_mappings[0];
         auto winner_mappings_pair = haplotype_ids.size() - 1 -
                 ordered_haplotype_barcodes_total_mappings[0];
         auto barcode_support_winners = std::make_pair(winner,
                                                       winner_pair);
         auto kmer_support_winners = std::make_pair(winner_mappings,
                                                    winner_mappings_pair);
         std::cout << "barcode winner: " << winner << " "
                   << winner_pair << std::endl;
         std::cout << "kmer winner: " <<winner_mappings<< " " << winner_mappings_pair
                   << std::endl;
         std::set<prm10xTag_t> s1;
         std::set<prm10xTag_t> s2;


         for (auto h:haplotype_ids[winner]) {
             std::cout << h << " ";
         }
             for (auto l: barcodes_supporting_haplotype[winner]) {
                 s1.insert(l);
             }
         for (auto l: barcodes_supporting_haplotype[winner_pair]) {
             s2.insert(l);
         }

         std::cout << std::endl;

         for (auto h:haplotype_ids[winner_mappings]) {
             std::cout << h << " ";
         }
         std::cout << std::endl;
         for (auto h:haplotype_ids[winner_mappings_pair]) {
             std::cout << oldnames[h] << " ";
         }
         std::cout << std::endl;
         for (auto h:haplotype_ids[winner_pair]) {
             std::cout << oldnames[h] << " ";
         }
         std::cout << std::endl;


         std::cout << std::endl;
         std::cout << "barcodes supporting size: " << s1.size() << " " << s2.size() << std::endl;
         if (s1.size() == 0){
             // these all seem to have the pattern aabbccddeeff......
             for (auto c: haplotype_barcodes_supporting){
                 std::cout << c << " ";
             }
             std::cout << std::endl;
             for (auto c: haplotype_barcodes_total_mappings){
                 std::cout << c << " ";
             }
             std::cout << std::endl;
         }
         // sometimes these are 0- why1?!?
         this->barcodes_supporting_winners = std::make_pair(s1, s2);
         std::cout << "barcodes supporting size: " << std::get<0>(this->barcodes_supporting_winners).size() << " "
                   << std::get<1>(this->barcodes_supporting_winners).size() << std::endl;

         if (winner== winner_mappings||
             winner == winner_mappings_pair) {
             if (winner_pair == winner_mappings ||
                     winner_pair == winner_mappings_pair) {
                 return 1;
             } else {
                 return 2;
             }
         } else {
             return 0;
         }
     }
     std::cout << " no relevant mappings" << std::endl;
     return 0;
 };

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



void HaplotypeScorer::find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >> bubbles){
    // algorithm: https://stackoverflow.com/questions/1867030/combinations-of-multiple-vectors-elements-without-repetition
    size_t P = 1;
    auto N =bubbles.size();
    for(size_t i=0;i<N;++i) {
        P *= bubbles[i].size();
        haplotype_nodes.insert(bubbles[i].begin(), bubbles[i].end());
    }
    std::vector<std::vector<sgNodeID_t >> haps;
    std::cout << P << " combinations to generate from " << N << " bubbles " << std::endl;
    for (size_t m=0; m < P; m++ ) {
        // this should hold the index to take from each bubble
        std::vector<size_t> indices(N);
        std::vector<sgNodeID_t > bubble;
        size_t m_curr = m;
        for (size_t i = 0; i < N; ++i) {
            indices[i] = m_curr% bubbles[i].size();
            bubble.push_back(bubbles[i][indices[i]]);
            m_curr /= bubbles[i].size();
        }

        haplotype_ids.push_back(bubble);
    }
    std::cout << haplotype_ids.size() << " haplotypes  generated " << std::endl;


    this->haplotype_nodes = haplotype_nodes;
    std::cout << "Haplotype nodes size: " << haplotype_nodes.size() << std::endl;

}