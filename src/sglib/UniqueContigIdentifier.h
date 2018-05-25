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
    std::vector<std::vector<sgNodeID_t > > GetAllUniqueContigs(double);
void WriteUniqueContentToFasta(int );
   void GetAllUniqueKmers();

private:

        SequenceGraph BuildGraph(std::string filename);
    bool UniqueKmersPerContig(int , sgNodeID_t,double );
    bool MatchingKmersPerContig(int , sgNodeID_t , std::vector<uint64_t > , double );
    bool CommonKmersPerContig(int , sgNodeID_t , double);

        int GetKmers(int);
    std::vector<SequenceGraph> graphs;
    std::map<uint64_t , std::set<int>> all_kmers;
    std::vector<std::vector<sgNodeID_t > > contig_ids;
    std::vector<std::string> assembly_names;
    std::vector<int>    kmers_per_assembly;
    uint64_t max_mem;

    std::vector<std::vector< uint64_t > > unique_kmers_for_assemblies;
    std::vector<std::vector< uint64_t > > common_kmers_for_assemblies;
};


#endif //BSG_UNIQUECONTIGIDENTIFIER_H
