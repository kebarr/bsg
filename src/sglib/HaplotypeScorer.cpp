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



/**
 *
 *
 * @brief
 * sum each barcode's score for each haplotype
 * @return
 * Returns a vector containing: total records generated by the factory,
 * total records after filtering,
 * total file records present in the file and the total number of records after filtering.
 */

void HaplotypeScorer::decide_barcode_haplotype_support(std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings, std::map<prm10xTag_t, std::vector<sgNodeID_t > > barcode_node_mappings){
    haplotype_barcodes_total_mappings.resize(haplotype_ids.size());
    haplotype_barcodes_supporting.resize(haplotype_ids.size());
    int total_mappings = 0;
    std::map<size_t , std::map< prm10xTag_t, int>> barcode_haplotype_shared;
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


    // previously found all nodes that a barcode maps to and only considered ones which mapped to enough of the haploype
    // for each haplotype, loop over each node and sum support
    for (int i = 0; i < haplotype_ids.size() ; i++){
        // to be same as previous, need to check each barcode counted maps to enough nodes in haplotype
        for (auto n: haplotype_ids[i]){
            if (node_tag_mappings[n].size() > 1) {

                for (auto f:node_tag_mappings[n]) {
                    // if this barcode shares
                    if (barcode_haplotype_shared[i][f.first] > node_tag_mappings[n][f.first]/2) {
                        haplotype_barcodes_supporting[i] += 1;

                            haplotype_barcodes_total_mappings[i] += f.second;
                            total_mappings += f.second;
                        barcodes_supporting_haplotype[i].push_back(f.first);

                    }
                }
            }
        }
    }


    std::cout << "Calculated haplotype support for each barcode, total mappings: " << total_mappings<<  std::endl;

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

     // need pairs with max barcodes + kmers support
     auto ordered_haplotype_barcodes_supporting = sort_indexes(haplotype_barcodes_supporting);
     for (auto o:ordered_haplotype_barcodes_supporting){
         std::cout << o << " " << haplotype_barcodes_supporting[o] << " ";
     }
     std::cout << std::endl;
     auto ordered_haplotype_barcodes_total_mappings = sort_indexes(haplotype_barcodes_total_mappings);
     for (auto o:ordered_haplotype_barcodes_total_mappings){
         std::cout << o << " " << haplotype_barcodes_total_mappings[o] << " ";
     }
     std::cout << std::endl;
     auto barcode_support_winners = std::make_pair(ordered_haplotype_barcodes_supporting[0], haplotype_ids.size() - 1 - ordered_haplotype_barcodes_supporting[0]);
     auto kmer_support_winners = std::make_pair(ordered_haplotype_barcodes_total_mappings[0], haplotype_ids.size() - 1 - ordered_haplotype_barcodes_total_mappings[0]);
    std::cout << "barcode winner: " << std::get<0>(barcode_support_winners) << " " << std::get<1>(barcode_support_winners) << std::endl;
     std::cout << "kmer winner: " << std::get<0>(kmer_support_winners) << " " << std::get<1>(kmer_support_winners) << std::endl;
     std::set<prm10xTag_t> s1;
     std::set<prm10xTag_t> s2;


     for (auto h:haplotype_ids[std::get<0>(barcode_support_winners)]){
        std::cout << h << " ";
         for (auto l: barcodes_supporting_haplotype[h]) {
             s1.insert(l);
         }
    }
     std::cout << std::endl;

     for (auto h:haplotype_ids[std::get<0>(kmer_support_winners)]){
         std::cout << h << " ";
     }
     std::cout << std::endl;
     for (auto h:haplotype_ids[std::get<1>(barcode_support_winners)]){
         //std::cout <<oldnames[h] << " ";
         for (auto l: barcodes_supporting_haplotype[h]) {
             s2.insert(l);
         }
     }
     std::cout << std::endl;

     for (auto h:haplotype_ids[std::get<1>(kmer_support_winners)]){
         std::cout << oldnames[h] << " ";
     }
     std::cout << std::endl;
     for (auto h:haplotype_ids[std::get<1>(barcode_support_winners)]){
         std::cout <<oldnames[h] << " ";
     }
     std::cout << std::endl;


     std::cout << std::endl;
     for (auto h:haplotype_ids[9]){
         std::cout <<oldnames[h] << " ";
     }
     std::cout << std::endl;

     for (auto h:haplotype_ids[6]){
         std::cout << oldnames[h] << " ";
     }
     std::cout << std::endl;
     barcodes_supporting_winners = std::make_pair(s1, s2);
     if (std::get<0>(barcode_support_winners) == std::get<0>(kmer_support_winners) || std::get<0>(barcode_support_winners) == std::get<1>(kmer_support_winners)){
         if ( std::get<1>(barcode_support_winners) == std::get<0>(kmer_support_winners) || std::get<1>(barcode_support_winners) == std::get<1>(kmer_support_winners) ){
             return 1;
         } else {
             return 2;
         }
     } else {
         return 0;
     }
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