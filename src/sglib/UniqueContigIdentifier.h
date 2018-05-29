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

struct Variant {
    SequenceGraph sg;
    std::set< uint64_t >  kmer_total;
    std::set< uint64_t >  unique_kmers;
    std::set< uint64_t >  common_kmers;
    std::set< uint64_t >  non_core_kmers;

    std::set<sgNodeID_t > common_contigs;
    std::vector<sgNodeID_t > unique_contigs;
    std::vector<sgNodeID_t > non_core_contigs;

    std::vector<double> shared_variation_all; // proportion of 'non common kmers' shared with each other var
    std::vector<double> shared_variation_non_core; // proportion of 'non common kmers' shared with each other var

};


class UniqueContigIdentifier {
public:
    void WriteAllUniqueContentToFasta();

        UniqueContigIdentifier(std::vector<std::string>, uint64_t );
    void FindCoreGenome(int);

   void GetUniqueContigs(int );
   void GetAllUniqueContigs(double);
    std::vector<std::vector<sgNodeID_t > > unique_contigs;
    void WriteUniqueContentToFasta(int , std::string);
   void GetAllUniqueKmers();
    void WriteAllNonCoreContentToFasta();
    void WriteNonCoreContentToFasta(int , std::string );



        void SummarizeVars();
private:
    void LoadGraphs();

        SequenceGraph BuildGraph(std::string filename);
    bool UniqueKmersPerContig(int , sgNodeID_t,double );
    bool MatchingKmersPerContig(int , sgNodeID_t , std::vector<uint64_t > , double );
    bool CommonKmersPerContig(int , sgNodeID_t , double);
    bool NonCoreKmersPerContig(int , sgNodeID_t , double match_thresh=0.5);

        size_t GetKmers(int);
    std::vector<Variant> variants;
    std::map<uint64_t , std::set<int>> all_kmers;
    std::vector<std::vector<sgNodeID_t > > contig_ids;
    std::vector<std::string> assembly_names;
    std::vector<int>    kmers_per_assembly;
    uint64_t max_mem;

    std::vector<std::vector< uint64_t > > unique_kmers_for_assemblies;
    std::vector<std::vector< uint64_t > > non_core_kmers_for_assemblies;

    std::vector< uint64_t >  common_kmers_for_assemblies;
};


#endif //BSG_UNIQUECONTIGIDENTIFIER_H
