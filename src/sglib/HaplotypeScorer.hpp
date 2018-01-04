//
// Created by Katie Barr (EI) on 14/11/2017.
//

#ifndef SG_HAPLOTYPE_SCORER_H
#define SG_HAPLOTYPE_SCORER_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <numeric>

#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>

class HaplotypeScorer{

public:
    // functions we will need:
    void find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >>);

    int score_haplotypes(std::vector<std::string>);
    std::vector<int > remove_nodes_with_no_barcode_support(std::vector< std::vector<prm10xTag_t> > tags_in_node, std::map<size_t , std::map< prm10xTag_t, int>>);
    int decide_barcode_haplotype_support(std::map<sgNodeID_t, std::map<prm10xTag_t, int > >, std::map<prm10xTag_t, std::vector<sgNodeID_t > >, std::vector< std::vector<prm10xTag_t> > );
    std::map<prm10xTag_t, std::map< int, int > > barcode_haplotype_mappings;
    bool success = false; // if doing partial success replace with enum
    std::pair<std::vector<sgNodeID_t >, std::vector< sgNodeID_t > > winners;
    //tags supporting h1 and h2
    std::pair<std::map<prm10xTag_t, int > , std::map<prm10xTag_t, int > > supporting_barcodes;
    std::pair<std::set<prm10xTag_t>, std::set<prm10xTag_t>  > barcodes_supporting_winners;


private:
    std::map<size_t , std::vector< prm10xTag_t> > barcodes_supporting_haplotype;
    // each het node
    std::set<sgNodeID_t > haplotype_nodes;
    // each possible hsplotype
    std::vector<std::vector<sgNodeID_t> > haplotype_ids;

    void print_voting_stats(std::vector<int> vote_totals);
    std::vector <prm10xTag_t> unused_barcodes;

    std::vector<int>  winner_for_barcode(prm10xTag_t barcode);
    // index is ha[plotype, value is number of barcodes supporting
    std::vector<int> haplotype_barcodes_supporting;
    // index is ha[plotype, value is number of kmers supporting

    std::vector<int> haplotype_barcodes_total_mappings;

};
#endif //SG_HAPLOTYPE_SCORER_H
