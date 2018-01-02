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

    void decide_barcode_haplotype_support(std::map<sgNodeID_t, std::map<prm10xTag_t, int > >, std::map<prm10xTag_t, std::vector<sgNodeID_t > >  );
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

    std::vector <prm10xTag_t> unused_barcodes;

    std::vector<int>  winner_for_barcode(prm10xTag_t barcode);
    std::vector<sgNodeID_t> haplotype_barcodes_supporting;
    std::vector<sgNodeID_t> haplotype_barcodes_total_mappings;

};
#endif //SG_HAPLOTYPE_SCORER_H
