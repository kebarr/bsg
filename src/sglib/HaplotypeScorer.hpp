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
    void load_haplotypes(std::string, int);

    void count_barcode_votes(PairedReadMapper &);
    int score_haplotypes(std::vector<std::string> );

    std::map<prm10xTag_t, std::map<sgNodeID_t , int> > barcode_node_mappings;
    void decide_barcode_haplotype_support();
    std::map<prm10xTag_t, std::map< int, int > > barcode_haplotype_mappings;
    bool success = false; // if doing partial success replace with enum
    std::pair<std::vector<sgNodeID_t >, std::vector< sgNodeID_t > > winners;
    //tags supporting h1 and h2
    std::pair<std::map<prm10xTag_t, int > , std::map<prm10xTag_t, int > > supporting_barcodes;

private:

    // each het node
    std::set<sgNodeID_t > haplotype_nodes;
    // each possible hsplotype
    std::vector<std::vector<sgNodeID_t> > haplotype_ids;

    std::vector <prm10xTag_t> unused_barcodes;

    std::vector<int>  winner_for_barcode(prm10xTag_t barcode);
    void analyse_scores(std::vector<std::string>, std::vector<int > , std::vector<int >, std::vector<int >, std::map<std::pair<int, int>, int> , std::map<std::pair<int, int>, int> , std::map<std::pair<int, int>, int>  );

    std::map<int, std::map<prm10xTag_t, int > > haplotype_barcode_agree;
    std::map<int, std::map<prm10xTag_t, int > > haplotype_barcode_disagree;

};
#endif //SG_HAPLOTYPE_SCORER_H
