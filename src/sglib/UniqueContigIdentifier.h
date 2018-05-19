//
// Created by Katie Barr (EI) on 19/05/2018.
//

#ifndef BSG_UNIQUECONTIGIDENTIFIER_H
#define BSG_UNIQUECONTIGIDENTIFIER_H
#include <string>
#include <map>
#include "sglib/SequenceGraph.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <sglib/KmerCompressionIndex.hpp>


class UniqueContigIdentifier {
public:
    UniqueContigIdentifier(std::vector<std::string>, uint64_t );
    void LoadGraphs();

    std::vector<sgNodeID_t > GetUniqueContigs(int );
    std::vector<sgNodeID_t > GetAllUniqueContigs(int );
void WriteUniqueContentToFasta(int );
    std::vector<sgNodeID_t > GetAllUniqueKmers();

private:
    SequenceGraph BuildGraph(std::string filename);
    void GetKmers(SequenceGraph &);
    std::vector<SequenceGraph> graphs;
    std::map<uint64_t , std::set<std::string>> all_kmers;
    std::vector<std::vector<sgNodeID_t > > contig_ids;
    std::map<std::string, int> assembly_names;
    uint64_t max_mem;
};


#endif //BSG_UNIQUECONTIGIDENTIFIER_H
