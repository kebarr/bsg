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

void HaplotypeScorer::decide_barcode_haplotype_support(){

    int support;
    int haplotypes_supported = 0;
    for (auto &mapping:barcode_node_mappings) {
        if (mapping.second.size() > 1) {
            std::vector<sgNodeID_t> nodes;
            std::vector<int> scores;
            for (auto n: mapping.second) {
                nodes.push_back(n.first);
                scores.push_back(n.second);
            }
            if (*std::max_element(scores.begin(), scores.end()) > 1) {
                for (int i = 0; i < haplotype_ids.size(); i++) {
                    std::vector<sgNodeID_t> nodes_in_haplotype;
                    std::vector<sgNodeID_t> h;
                    h = haplotype_ids[i];
                    // find all nods in each haplotype that this barcode maps to
                    for (auto n1: nodes) {
                        if (std::find(h.begin(), h.end(), n1) != h.end()) {
                            nodes_in_haplotype.push_back(n1);
                        }
                    }
                    // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                    if (nodes_in_haplotype.size() >= (nodes.size() / 2) && nodes_in_haplotype.size() > 1) {
                        support = 0;
                        for (auto a: nodes_in_haplotype) {
                            support += mapping.second[a];
                        }
                        barcode_haplotype_mappings[mapping.first][i] = support;
                        haplotypes_supported += 1;
                    } else {
                        unused_barcodes.push_back(mapping.first);
                    }
                }
            }

        } else {
            unused_barcodes.push_back(mapping.first);
        }
        haplotypes_supported = 0;

    }

    std::cout << "Calculated haplotype support for each barcode, " << barcode_haplotype_mappings.size() <<  std::endl;

}


/**
 *
 *
 * @brief
 * Map reads to graph, sum up matches for each comntig mapped to by each barcode  barcode
 * @return
 * Load read files 1 and 2, populating map of id to barcode (this should be a vector the same as read to node but either my input reads are dodgy or there's something else wrong)
 * find mappings for each node in the bubble contigs. for each barcode, sum up the number of mappings for each node that a read with that barcode maps to
 *
 */
void HaplotypeScorer::count_barcode_votes(PairedReadMapper & mapper){
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << " nodes" << std::endl;
    std::cout << "NOw counting barcode votes... " << std::endl;
    int counter = 0;
    for (auto r: mapper.reads_in_node){
        counter += 1;
    }
    for (auto &node:haplotype_nodes){
        for (auto &mapping:mapper.reads_in_node[node>0?node:-node]){
            auto barcode = mapper.read_to_tag[mapping.read_id];
             barcode_node_mappings[barcode][node] += mapping.unique_matches;
            counter += 1;
        }
    }
    std::cout << "counted " << counter << " votes for " << barcode_node_mappings.size() << "barcodes" << std::endl;
};

/**
 *
 *
 * @brief
 * Decide winning phasing
 * @return
 * Calculate, for each haplotype, a number of metrics that represent how well the 10x reads support that phasing
 * loop over each barcode in turn, then each phasing pair
 * \todo decide criteria for winning haploype
 *
 */
int HaplotypeScorer::score_haplotypes(std::vector<std::string> oldnames) {
    /*std::cout << "{{";
    for (auto h: haplotype_ids){
        for (auto n: h){
            std::cout<<"\"" << oldnames[n] << "\", ";
        }
        std::cout<< "}, \n {";
    }
    std::cout << "}}\n";*/

    auto number_haplotypes = haplotype_ids.size();

    std::cout << "Finding most supported of " << number_haplotypes<< " possible haplotypes"<<std::endl;
    //initialize score arrays- index is haplotype index
    std::vector<int > haplotype_support;
    std::vector<int > haplotype_not_support;
    std::vector<int> haplotype_overall_support;
    haplotype_support.resize(number_haplotypes);
    haplotype_not_support.resize(number_haplotypes);
    haplotype_overall_support.resize(number_haplotypes);
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int> hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;
    prm10xTag_t barcode;
    for (auto &bm: barcode_haplotype_mappings) {
        barcode = bm.first;
        std::vector<int> winners = winner_for_barcode(barcode); // ideally should be length 1
        for (auto winner:winners){
            int pair = haplotype_ids.size() - 1 - winner;
            std::pair<int, int> res = pair > winner ? std::make_pair(winner, pair) : std::make_pair( pair, winner);
            // like this, support and pair support always identical
            haplotype_support[winner] += 1;
            hap_pair_support[res] += 1;
            haplotype_barcode_agree[winner][barcode] += bm.second[winner];
            haplotype_barcode_disagree[winner][barcode] += bm.second[pair];
        }
        for (int hap = 0; hap < number_haplotypes/ 2; hap++) {
            // pair = len(self.list_of_possible_haplotypes) -1 -haplotype
            auto pair = number_haplotypes - 1 - hap;
            if (bm.second.find(hap) != bm.second.end()) {
                haplotype_overall_support[hap] += bm.second[hap];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[hap];
            }

            if (bm.second.find(pair) != bm.second.end()) {
                haplotype_overall_support[pair] += bm.second[pair];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[pair];
            }
            if (bm.second.find(hap) == bm.second.end()) {
                haplotype_not_support[hap] += 1;
            }
            if (bm.second.find(pair) == bm.second.end()) {
                haplotype_not_support[pair] += 1;
            }
            if (bm.second.find(hap) == bm.second.end() and bm.second.find(pair) == bm.second.end()) {
                hap_pair_not_support[std::make_pair(hap, pair)] += 1;

            }
        }
    }
    //analyse_scores(oldnames, haplotype_support, haplotype_not_support, haplotype_overall_support, hap_pair_support,hap_pair_not_support, hap_pair_support_total_score );

    auto support_max_index = std::distance(haplotype_support.begin(), std::max_element (haplotype_support.begin(),haplotype_support.end()));
    auto overall_support_max_index = std::distance(haplotype_overall_support.begin(), std::max_element (haplotype_overall_support.begin(),haplotype_overall_support.end()));
    // these just give last....
    //auto pair_support_max_index = std::max_element (std::begin(hap_pair_support),std::end(hap_pair_support));
    //auto pair_overall_support_max_index = std::max_element (std::begin(hap_pair_support_total_score),std::end(hap_pair_support_total_score));
    int pair_support_max = 0;
    std::pair<sgNodeID_t , sgNodeID_t > pair_support_winner;
    int pair_overall_support_max = 0;
    std::pair<sgNodeID_t , sgNodeID_t > pair_support_overall_winner;

    for (auto p: hap_pair_support){
        if (p.second > pair_support_max){
            pair_support_max = p.second;
            pair_support_winner = std::make_pair(std::get<0>(p.first), std::get<1>(p.first));
        }
    }
    for (auto p: hap_pair_support_total_score){
        if (p.second > pair_overall_support_max){
            pair_overall_support_max = p.second;
            pair_support_overall_winner = std::make_pair(std::get<0>(p.first), std::get<1>(p.first));
        }
    }
    std::cout << "Support max index: " << support_max_index << " max support value: " << haplotype_support[support_max_index] << std::endl;
    std::cout << "overall Support max index: " << overall_support_max_index << " max support value: " << haplotype_overall_support[overall_support_max_index] << std::endl;
    std::cout << "pair Support max index: " << std::get<0>(pair_support_winner) << " " << std::get<1>(pair_support_winner)<< " " << " max support value: " << pair_support_max << std::endl;
    std::cout << "pair overall Support max index: " << std::get<0>(pair_support_overall_winner) << " " << std::get<1>(pair_support_overall_winner)<< " " << " max support value: " << pair_overall_support_max   << std::endl;

    std::cout << "Haplotype support:\n";
    print_int_vector(haplotype_support);

    std::cout << "Haplotype overall support:\n";
    print_int_vector(haplotype_overall_support);

    std::cout << "Pair support \n";
    print_pair_int_map(hap_pair_support);


    std::cout << "Pair overall support \n";
    print_pair_int_map(hap_pair_support_total_score);
    // stop now... but TODO: sum supports and overall supports, see if they vary each time
    //saw no variations when loading, only dumping
    // also see if it varies if comment out parallel code- seems not to when dumping
    std::cout << "std::get<0>(pair_support_overall_winner) " << std::get<0>(pair_support_overall_winner) << "  haplotype_overall_support[overall_support_max_index]  " <<  haplotype_overall_support[overall_support_max_index] << " overall_support_max_index " << overall_support_max_index<< " std::get<1>(pair_support_overall_winner) "<< std::get<1>(pair_support_overall_winner) << "  std::get<0>(pair_support_winner) " << std::get<0>(pair_support_winner) << " std::get<1>(pair_support_winner) " << std::get<1>(pair_support_winner) << " haplotype_support[support_max_index]) " << haplotype_support[support_max_index] << " support_max_index " << support_max_index << std::endl;

    if ((std::get<0>(pair_support_overall_winner) == overall_support_max_index ||
         std::get<1>(pair_support_overall_winner) == overall_support_max_index ) &&
        (std::get<0>(pair_support_winner) == support_max_index ||
         std::get<1>(pair_support_winner) == support_max_index )) {
        this->winners = std::make_pair(haplotype_ids[std::get<0>(pair_support_overall_winner)], haplotype_ids[std::get<1>(pair_support_overall_winner)]);
        this->success = true;
        return  1;
    }

    return 0;
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

void HaplotypeScorer::analyse_scores(std::vector<std::string > oldnames, std::vector<int > haplotype_support, std::vector<int > haplotype_not_support, std::vector<int > haplotype_overall_support, std::map<std::pair<int, int>, int> hap_pair_support, std::map<std::pair<int, int>, int>  hap_pair_not_support, std::map<std::pair<int, int>, int> hap_pair_support_total_score ) {
    std::vector<int> haplotype_support_vals;
    std::vector<int> haplotype_not_support_vals;
    std::vector<int> haplotype_overall_support_vals;
    for (int i = 0; i < haplotype_ids.size(); i++) {
        haplotype_support_vals.push_back(haplotype_support[i]);
        haplotype_not_support_vals.push_back(haplotype_not_support[i]);
        haplotype_overall_support_vals.push_back(haplotype_overall_support[i]);

    }

    std::cout << "overall support vals: \n";

    print_int_vector(haplotype_overall_support_vals);
    std::cout << "support vals: \n";

    print_int_vector(haplotype_support_vals);
    // success is judged by getting highest scoring pair/non support/overall winners and checking they agree
    auto support_max_index = std::distance(haplotype_support.begin(), std::max_element (haplotype_support.begin(),haplotype_support.end()));
    auto overall_support_max_index = std::distance(haplotype_overall_support.begin(), std::max_element (haplotype_overall_support.begin(),haplotype_overall_support.end()));
    // these just give last....
    //auto pair_support_max_index = std::max_element (std::begin(hap_pair_support),std::end(hap_pair_support));
    //auto pair_overall_support_max_index = std::max_element (std::begin(hap_pair_support_total_score),std::end(hap_pair_support_total_score));
    int pair_support_max = 0;
    std::pair<sgNodeID_t , sgNodeID_t > pair_support_winner;
    int pair_overall_support_max = 0;
    std::pair<sgNodeID_t , sgNodeID_t > pair_support_overall_winner;

    for (auto p: hap_pair_support){
        if (p.second > pair_support_max){
            pair_support_max = p.second;
            pair_support_winner = std::make_pair(std::get<0>(p.first), std::get<1>(p.first));
        }
    }
    for (auto p: hap_pair_support_total_score){
        if (p.second > pair_overall_support_max){
            pair_overall_support_max = p.second;
             pair_support_overall_winner = std::make_pair(std::get<0>(p.first), std::get<1>(p.first));
        }
    }
    std::cout << "Support max index: " << support_max_index << " max support value: " << haplotype_support[support_max_index] << std::endl;
    std::cout << "overall Support max index: " << overall_support_max_index << " max support value: " << haplotype_overall_support[overall_support_max_index] << std::endl;
    std::cout << "pair Support max index: " << std::get<0>(pair_support_winner) << " " << std::get<1>(pair_support_winner)<< " " << " max support value: " << pair_support_max << std::endl;
    std::cout << "pair overall Support max index: " << std::get<0>(pair_support_overall_winner) << " " << std::get<1>(pair_support_overall_winner)<< " " << " max support value: " << pair_overall_support_max   << std::endl;

    std::cout << "NNNNNN*****&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n";
    std::vector<std::pair<int, int> > supports;
    std::vector<std::pair<int, int> > overall_supports;
    std::vector<std::pair<std::pair<int, int>, int> > pair_supports;
    std::vector<std::pair<std::pair<int, int>, int> > pair_overall_supports;
    for (int i = 0; i < haplotype_ids.size(); i++) {
        supports.push_back(std::make_pair(i, haplotype_support[i]));
        overall_supports.push_back(std::make_pair(i, haplotype_support[i]));
    }
    std::sort(supports.begin(), supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    std::sort(overall_supports.begin(), overall_supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    // when not rushing, get rid of al this repetition
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(),
                                                haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);


    std::cout << "Haplotype support max : " << *support_max << " min: "
              << *std::min_element(haplotype_support_vals.begin(), haplotype_support_vals.end()) << " mean: "
              << support_mean << " stdev: " << support_stdev << std::endl;
    std::cout << "Haplotype not support max : "
              << *std::max_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " min: "
              << *std::min_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " mean: "
              << avg(haplotype_not_support_vals) << std::endl;
    std::cout << "Haplotype overall support max: " << *overall_support_max << " min: "
              << *std::min_element(haplotype_overall_support_vals.begin(), haplotype_overall_support_vals.end())
              << " mean: " << overall_support_mean << " stdev: " << overall_stdev << std::endl;
    if (hap_pair_support.size() > 0 && hap_pair_support_total_score.size() > 0) {
        std::vector<int> hap_pair_not_support_values;
        std::vector<int> hap_pair_support_values;
        for (auto h : hap_pair_support) {
            hap_pair_support_values.push_back(h.second);
            pair_supports.push_back(std::make_pair(h.first, h.second));
        }
        std::cout << "overall support vals: \n";

        print_int_vector(haplotype_overall_support_vals);
        std::cout << "support vals: \n";

        print_int_vector(haplotype_support_vals);
        std::cout << "pair support vals: \n";

        print_int_vector(hap_pair_support_values);

        for (auto h : hap_pair_not_support) {
            hap_pair_not_support_values.push_back(h.second);
        }
        std::vector<int> hap_pair_support_total_score_values;
        for (auto h : hap_pair_support_total_score) {
            hap_pair_support_total_score_values.push_back(h.second);
            pair_overall_supports.push_back(std::make_pair(h.first, h.second));
        }

        std::cout << "pair support vals: \n";

        print_int_vector(hap_pair_support_total_score_values);
        std::sort(pair_supports.begin(), pair_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });
        std::sort(pair_overall_supports.begin(), pair_overall_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });

        auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());
        auto not_pair_support_max = std::max_element(hap_pair_not_support_values.begin(),
                                                     hap_pair_not_support_values.end());
        auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(),
                                                         hap_pair_support_total_score_values.end());
        auto pair_support_min = std::min_element(hap_pair_support_values.begin(), hap_pair_support_values.end());
        auto not_pair_support_min = std::min_element(hap_pair_not_support_values.begin(),
                                                     hap_pair_not_support_values.end());
        auto overall_pair_support_min = std::min_element(hap_pair_support_total_score_values.begin(),
                                                         hap_pair_support_total_score_values.end());
        auto pair_support_mean = avg(hap_pair_support_values);
        auto pair_support_stdev = stdev(hap_pair_support_values, pair_support_mean);
        auto pair_not_support_mean = avg(hap_pair_not_support_values);
        auto pair_not_support_stdev = stdev(hap_pair_not_support_values, pair_not_support_mean);
        auto pair_overall_support_mean = avg(hap_pair_support_total_score_values);
        auto pair_overall_support_stdev = stdev(hap_pair_support_total_score_values, pair_overall_support_mean);
        std::cout << "pair support size: " << pair_supports.size() << " total: "
                  << hap_pair_support_total_score_values.size() << std::endl;
        std::cout << "Haplotype pair support max : " << *pair_support_max << "min : " << *pair_support_min << " mean: "
                  << pair_support_mean << " stdev: " << pair_support_stdev << std::endl;
        std::cout << "Haplotype pair not support max : " << *not_pair_support_max << "min : " << *not_pair_support_min
                  << " mean: " << pair_not_support_mean << " stdev: " << pair_not_support_stdev << std::endl;
        std::cout << "Haplotype pair overall support max : " << *overall_pair_support_max << "min : "
                  << *overall_pair_support_min << " mean: " << pair_overall_support_mean << " stdev: "
                  << pair_overall_support_stdev << std::endl;
        // get winners
        std::vector<int> support_winner;
        std::vector<int> overall_support_winner;
        for (int h = 0; h < haplotype_ids.size(); h++) {
            if (haplotype_support[h] == *support_max) {
                support_winner.push_back(h);
            }
            if (haplotype_overall_support[h] == *overall_support_max) {
                overall_support_winner.push_back(h);
            }
        }
        std::vector<std::pair<int, int> > pair_support_winner;
        std::vector<std::pair<int, int> > pair_overall_support_winner;

        for (auto h: hap_pair_support) {
            if (h.second == *pair_support_max) {
                pair_support_winner.push_back(h.first);
            }
        }
        for (auto h: hap_pair_support_total_score) {
            if (h.second == *overall_pair_support_max) {
                pair_overall_support_winner.push_back(h.first);

            }
        }
        std::cout << "Support winner: ";
        print_int_vector(support_winner);
        std::cout << "overall SUpport winner: ";
        print_int_vector(overall_support_winner);
        std::cout << "pair SUpport winner: ";
        print_pair_int_vector(pair_support_winner);
        std::cout << "pair overall SUpport winner: ";
        for (auto w:pair_overall_support_winner){
            std::cout << oldnames[std::get<0>(w)] << " " << oldnames[std::get<1>(w)] << " ";
        }
        std::cout << std::endl;
        print_pair_int_vector(pair_overall_support_winner);
        for (auto a: haplotype_ids[support_winner[0]]){
            std::cout << oldnames[a]     << " ";
        }
        std::cout << std::endl;
        print_id_vector(haplotype_ids[support_winner[0]]);


    }
}

/**
 *
 *
 * @brief
 * Find haplotype that a barcode votes for
 * * @return
 * loop over barcode aplotype mappings to find maximum support
 * if more than2 have same support both are returned
 * \todo decide heuristics to vote for winner
 */
std::vector<int>  HaplotypeScorer::winner_for_barcode(prm10xTag_t barcode){
    int max=0;
    std::vector<int> winners;
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second > max){
            max = h.second;
        }
    }
    //TODO: DECIDE CRITERIA FOR MINIMUM SUPPORT
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second == max){
            winners.push_back(h.first);
        }
    }
    return winners;
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