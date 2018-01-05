//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "HaplotypeScorer.hpp"

/*
 * find bubbles with mappings that meet minimum phaseability requirements
 * calculate possible haplotypes
 * find mapping support for each haplotype
 * a barcode supports a haplotype h1 if it maps to at least 2 nodes in h1 and twice as many nodes from h1 as it does h2
 * for each haplotype, sum barcode support = kmer support
 *
 *
 */

HaplotypeScorer::HaplotypeScorer(std::vector<sgNodeID_t> component) : component(component){};


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


std::vector<std::vector<sgNodeID_t >> HaplotypeScorer::find_supported_nodes(std::vector<std::vector<sgNodeID_t >> bubbles, std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings, std::map<prm10xTag_t, std::set<sgNodeID_t > > barcode_node_mappings){
    // should find supported nodes before calculating haplotypes
    // bubbles already filtered so only ones which have at least 1 contig with barcode mappings included
    std::vector<sgNodeID_t > phaseable;
    std::vector<sgNodeID_t > not_phaseable;
    // bubbles where max 1 node isn't phaseable should be phaseable
    std::vector<std::vector<sgNodeID_t >> phaseable_bubbles;
    std::cout << "Finding nodes phased by mapper " << std::endl;
    // whether node is supported is haplotype independent- node can be phased by 10x if reads map to that and another node
    for (auto bubble:bubbles){
        int phaseable_for_bubble = 0;
        for (auto node: bubble) {
            std::cout << node << ": ";
            // thias  should always contain at least 2 mappnig
            auto tags_mapped_to = node_tag_mappings[node];
            bool confirmed_phaseable = false;
            for (auto t:tags_mapped_to) {
                // majority of barcodes map to single node
                std::cout << t.first << " maps to node:  " << t.second << "times, and maps to other nodes:";
                for (auto b: barcode_node_mappings[t.first]){

                    std::cout << " " << b << " ";
                }
                std::cout << std::endl;

                if (barcode_node_mappings[t.first].size() > 1) {
                    std::cout << "phaseable " << std::endl;
                    // node is mapped to barcode which maps to other nodes, so is phaseable
                    phaseable.push_back(node);
                    phaseable_for_bubble += 1;
                    // know this node is phaseable, so ove to next
                    confirmed_phaseable = true;
                    break;

                } else {
                    not_phaseable.push_back(node);
                }
            }
            if (confirmed_phaseable){
                break;
            }
        }
        if (phaseable_for_bubble >= bubble.size() -1){
            phaseable_bubbles.push_back(bubble);
        }
    }
    std::cout << phaseable.size() << " nodes from " << haplotype_nodes.size() << " phaseable by mapping " << std::endl;
    std::cout << phaseable_bubbles.size() << " bubbles from " << bubbles.size() << " phaseable by mapping " << std::endl;
    return  phaseable_bubbles;
}


std::map<size_t , std::map< prm10xTag_t, int>> HaplotypeScorer::count_barcode_haplotype_mappings(std::map<prm10xTag_t, std::set<sgNodeID_t > > barcode_node_mappings){
    std::map<size_t , std::map< prm10xTag_t, int>> barcode_haplotype_shared;// number of nodes that a barcode maps to in each hap

    // for each barcode, calculate which haps it shares nodes with
    for (auto m: barcode_node_mappings){
        auto barcode = m.first;
        auto nodes = m.second;
        for (int i = 0; i < haplotype_ids.size() ; i++) {
            for (auto l: nodes) {
                if (std::find(haplotype_ids[i].begin(), haplotype_ids[i].end(), l) != haplotype_ids[i].end()) {
                    // this is only ever 1 in my stdout, but should equal number of nodes in haplotype that barcode maps to
                    barcode_haplotype_shared[i][barcode] += 1;

                }

            }
        }
    }
    return  barcode_haplotype_shared;
};

/**
 *
 *
 * @brief
 * sum each barcode's score for each haplotype
 * @return
 * Returns number of barcodes that can be used to phase this component- i.e. map to more than 2 het sites
 * For each barcode, calculate how many nodes it shares with each haplotype- i.e. how many nodes it ,aps to contained in each haplotype
 * then for each node in each haplotype, find barcodes mapping to that node which do not map to alternate haplotype nodes
 */


int HaplotypeScorer::decide_barcode_haplotype_support(std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings_in, std::map<prm10xTag_t, std::set<sgNodeID_t > > barcode_node_mappings){
    std::cout << "calculating support from each barcodes for " << haplotype_ids.size() << " haplotypes \n";
    node_tag_mappings =node_tag_mappings_in;
    //auto to_ignore = remove_nodes_with_no_barcode_support(node_tag_mappings, barcode_haplotype_shared);
    // updating haplotypes changes pair indexing!!! aso don't, make list of haplotypes to skip
    haplotype_barcodes_total_mappings.resize(haplotype_ids.size());
    haplotype_barcodes_supporting.resize(haplotype_ids.size());
    int total_mappings = 0;
    int total_barcodes = 0;
    auto barcode_haplotype_shared = count_barcode_haplotype_mappings(barcode_node_mappings);// number of nodes that a barcode maps to in each hap
    size_t shared_nodes = 0;




    // previously found all nodes that a barcode maps to and only considered ones which mapped to enough of the haploype
    // for each haplotype, loop over each node and sum support
    for (int i = 0; i < haplotype_ids.size() ; i++){
            int pair = haplotype_ids.size() - i - 1;
            // to be same as previous, need to check each barcode counted maps to enough nodes in haplotype
                    for (auto barcode_tot : barcode_haplotype_shared[i]) {
                        auto barcode = barcode_tot.first;
                        //sum barcodes supporting each haplotype
                        if (barcode_haplotype_shared[i][barcode] > barcode_haplotype_shared[pair][barcode] * 2) {
                            haplotype_barcodes_supporting[i] += 1;
                            for (auto node: haplotype_ids[i]) {
                                // for each node i hap, add mappings from reads with this barcode
                                //if (node_tag_mappings[node][barcode] > ) {// if more than 1 barcode maps tp this node
                                    haplotype_barcodes_total_mappings[i] += 1;// node_tag_mappings[node][barcode];
                                    // number of kmers mapped to that node[e from reads with that barcode
                                    total_mappings += 1;//node_tag_mappings[node][barcode];
                                //std::cout << "node_tag_mappings[node][barcode];" << node << " " << barcode << " " << node_tag_mappings[node][barcode]<< "\n";
                                //}
                            }
                                std::cout << " nnnnnnnnnnnnnnnn \n";
                                std::cout << "barcode contributing: " << barcode << " " << barcode_node_mappings[barcode].size() << " to hap " << i << std::endl;
                                for (auto k:barcode_node_mappings[barcode]){
                                    std::cout << k << " ";
                                }
                                std::cout << "\n";
                            total_barcodes += 1;
                            if (barcodes_supporting_haplotype[i].size() > 0) {
                                barcodes_supporting_haplotype[i].push_back(barcode);
                            } else {
                                barcodes_supporting_haplotype[i] = {barcode};

                            }

                    }

                }
            }

        //}
for (int i=0; i < haplotype_barcodes_supporting.size() ; i++){
    if (haplotype_barcodes_supporting[i] > 0 || haplotype_barcodes_total_mappings[i] > 0) {
        //auto ignored = std::find(to_ignore.begin(), to_ignore.end(), i) == to_ignore.end();
        std::cout << "i: " << i << " votes " << haplotype_barcodes_supporting[i] << " kmers "
                  << haplotype_barcodes_total_mappings[i] << std::endl;
        int sum = 0;
        for (auto node: haplotype_ids[i]) {
            sum += node_tag_mappings[node].size();
        }
        std::cout << std::endl << " sum: " << sum << std::endl;
    }
}
    std::cout << haplotype_barcodes_supporting[0] << " " << haplotype_barcodes_supporting[1] << std::endl;
    std::cout << haplotype_barcodes_total_mappings[0] << " " << haplotype_barcodes_total_mappings[1] << std::endl;


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

int HaplotypeScorer::decide_results(std::vector<std::string> oldnames, std::pair< sgNodeID_t, sgNodeID_t >  barcode_support_winners, std::pair< sgNodeID_t, sgNodeID_t >kmer_support_winners){

    auto winner = std::get<0>(barcode_support_winners);
    auto winner_pair = std::get<1>(barcode_support_winners);

    auto winner_mappings = std::get<0>(kmer_support_winners);
    auto winner_mappings_pair = std::get<1>(kmer_support_winners);

    std::cout << "barcode winner: " << winner << " "
              << winner_pair << std::endl;
    std::cout << "kmer winner: " << winner_mappings << " " << winner_mappings_pair
              << std::endl;
    std::set <prm10xTag_t> s1;
    std::set <prm10xTag_t> s2;


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
    if (s1.size() == 0) {
        for (auto c: haplotype_barcodes_supporting) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        for (auto c: haplotype_barcodes_total_mappings) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    barcodes_supporting_winners = std::make_pair(s1, s2);
    std::cout << "barcodes supporting size: " << std::get<0>(barcodes_supporting_winners).size() << " "
              << std::get<1>(barcodes_supporting_winners).size() << std::endl;

    if (winner == winner_mappings ||
        winner == winner_mappings_pair) {
        if (winner_pair == winner_mappings ||
            winner_pair == winner_mappings_pair) {
            return 1;
        } else {
            return 2;
        }
    } else if (winner == winner_mappings ||
               winner == winner_mappings_pair){
        return  2;
    } else {
        return 0;
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


    // if there was a vote (rather than 0 barcodes for each haplotype) and no tie, check winner
     if (haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[0]] >0 && haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[0]] > 0) {

         if (haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[0]] !=
             haplotype_barcodes_supporting[ordered_haplotype_barcodes_supporting[1]] &&
             haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[0]] !=
             haplotype_barcodes_total_mappings[ordered_haplotype_barcodes_total_mappings[1]]
                 ) {
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
             return decide_results(oldnames, barcode_support_winners, kmer_support_winners);

         } else {

             auto winners_supporting = select_nodes_from_draw(ordered_haplotype_barcodes_supporting);
             if (std::get<0>(winners_supporting)[0] == 0){
                 return 0;
             }
             auto winners_total_supporting = select_nodes_from_draw(ordered_haplotype_barcodes_total_mappings);
             if (std::get<0>(winners_total_supporting)[0] == 0){
                 return 0;
             }
             size_t match_len = 0;
             size_t  par_match_len =0;
             int limit = std::get<0>(winners_supporting).size() >= std::get<0>(winners_total_supporting).size() ? std::get<0>(winners_supporting).size(): std::get<0>(winners_total_supporting).size();
             for (int w=0; w < limit; w++){
                 auto w0 = std::get<0>(winners_supporting)[w];
                 auto w1 = std::get<1>(winners_supporting)[w];
                 auto t0 = std::get<0>(winners_total_supporting)[w];
                 auto t1 = std::get<1>(winners_total_supporting)[w];

                 std::set <prm10xTag_t> s1;
                 std::set <prm10xTag_t> s2;
                 for (auto l: node_tag_mappings[w0]) {
                     s1.insert(l.first);
                 }
                 for (auto l: node_tag_mappings[w1]) {
                     s2.insert(l.first);
                 }
                 barcodes_supporting_winners = std::make_pair(s1, s2);

             if (w0 == t0 && w1 == t1 ) {
                 if (w1==t0 && w0 == t1){

                     match_len += 1;
                 } else {
                     par_match_len += 1;
                 }

                } else if (w1==t0 && w0 == t1){
                 par_match_len += 1;
                }
             }
             if (match_len == limit ) {
                 return 1;
             } else if (par_match_len == limit){
                 return  2;
             } else {
                 return  0;
             }


         }
     }
     std::cout << " no relevant mappings" << std::endl;
     return 0;
 };


std::pair < std::vector<sgNodeID_t >, std::vector<sgNodeID_t > > HaplotypeScorer::select_nodes_from_draw(std::vector<size_t > order_winners){
   // can potentially get soomething if draw...
    // if 2 nodes fro same bubble, can't

    std::cout << order_winners[0] << " won with score " << haplotype_barcodes_supporting[order_winners[0]] << " " << order_winners[1] << " runner up with score " << haplotype_barcodes_supporting[order_winners[1]] << std::endl;
for (auto h: haplotype_ids[order_winners[0]]){
    std::cout << h << " ";
}
    std::cout << std::endl;
    //if its a draw, take all contigs in common from top equals
    int top_score = haplotype_barcodes_supporting[order_winners[0]];
    std::vector<size_t > top_phasings;
    std::set<size_t > bubbles_seen;
    top_phasings.push_back(order_winners[0]);
    for (size_t s=1; s < order_winners.size(); s++){
        auto score = haplotype_barcodes_supporting[order_winners[s]];
        if (score==top_score){
            top_phasings.push_back(s);
            std::cout << "s: "<< s << " os: " <<order_winners[s] << " score: " << score << std::endl;
            for (auto n:haplotype_ids[s]){
                auto bubble = bubble_map[n];
                auto size1 = bubbles_seen.size();
                bubbles_seen.insert(bubble);
                auto size2 = bubbles_seen.size();
                if (size1 == size2){
                    std::cout << "ambiguous draw \n";
                    sgNodeID_t r = 0;
                    std::vector<sgNodeID_t > res;
                    res.push_back(r);
                    return  std::make_pair(res, res);
                }

            }
        } else {
            break;
        }
    }
    std::cout << top_phasings.size() << " drawn top phasings"<< std::endl;
    std::vector<sgNodeID_t > grouped_nodes;
    std::set<sgNodeID_t > grouped_nodes_set;
    // need to ensure that we don't take two from same bubble
    for (int p=0; p < top_phasings.size()-1; p++) {
        for (auto node:haplotype_ids[p]){
            std::cout << "node " << node << " m: ";
            for (auto m: node_tag_mappings[node]) {
                std::cout << m.first << " " << m.second << " ";
            }
            std::cout << std::endl;
        }
        // this assumes that both haps don't draw
        std::set_intersection(haplotype_ids[p].begin(), haplotype_ids[p].end(), haplotype_ids[p+1].begin(), haplotype_ids[p+1].end(),std::inserter(grouped_nodes, grouped_nodes.begin()));
        std::cout << "p:  << " << p << "intersected nodes: " << grouped_nodes.size() << " ";
        for (auto n:grouped_nodes){
            std::cout << n << " ";
        }
        std::cout << std::endl;
    }
    //TODO:don't assume diploid!!
    std::vector<sgNodeID_t > pair_grouped_nodes;
    for (auto w:grouped_nodes){
        auto bubble = bubbles[bubble_map[w]];
        for (auto node: bubble) {
            if (node != w) {
                pair_grouped_nodes.push_back(node);
                //std::cout << " node: " << w << " pair: " << node ;
                break;
            }
        }
        std::cout << std::endl;

    }

    return std::make_pair(grouped_nodes, pair_grouped_nodes);
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



void HaplotypeScorer::find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >> bubbles_in, std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings, std::map<prm10xTag_t, std::set<sgNodeID_t > > barcode_node_mappings){
    // only use bubbles which meet miinimum requirements to be phaseable- nodes mapped to tag which is mapped to at least 1 other node
    auto bubbles_supported = find_supported_nodes(bubbles_in, node_tag_mappings, barcode_node_mappings);
    // above results in nothing being phaseable- because lots of barcodes only map to single node
    //auto bubbles_supported = bubbles_in;
    // fr time being so that i can use these reads for next bit
    bubbles_supported = bubbles_in;
    bubbles = bubbles_in;
    if (bubbles_supported.size() > 1) {
        // algorithm: https://stackoverflow.com/questions/1867030/combinations-of-multiple-vectors-elements-without-repetition
        size_t P = 1;
        auto N = bubbles_supported.size();
        for (size_t i = 0; i < N; ++i) {
            P *= bubbles_supported[i].size();
            haplotype_nodes.insert(bubbles_supported[i].begin(), bubbles_supported[i].end());
            for (auto n: bubbles_supported[i]) {
                bubble_map[n] = i;
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
                indices[i] = m_curr % bubbles_supported[i].size();
                bubble.push_back(bubbles_supported[i][indices[i]]);
                m_curr /= bubbles_supported[i].size();
            }

            haplotype_ids.push_back(bubble);
        }
        std::cout << haplotype_ids.size() << " haplotypes  generated " << std::endl;

        for (auto n:component){
            if (std::find(haplotype_nodes.begin(), haplotype_nodes.end(), n) == haplotype_nodes.end()){
                hom_nodes.push_back(n);
            }
        }
        std::cout << "Haplotype nodes size: " << haplotype_nodes.size() << std::endl;
    }


}