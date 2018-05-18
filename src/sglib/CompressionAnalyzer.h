//
// Created by Katie Barr (EI) on 02/03/2018.
//

#ifndef BSG_COMPRESSIONANALYZER_H
#define BSG_COMPRESSIONANALYZER_H

#include <string>
#include <map>
#include "sglib/SequenceGraph.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/CoreGenomeFinder.h>
#include <sglib/RepeatAnalyzer.h>


struct NodeCompressions {
    std::string lib_name_r1;
    std::string lib_name_r2;
    int index;

    //std::map<sgNodeID_t, double> compressions; - just use sg indexing
    std::vector< double> compressions;
    std::vector<std::vector<sgNodeID_t >> canonical_repeats;
};

/*
struct Repeat{
    sgNodeID_t repeat;
    std::vector<sgNodeID_t> in;
    std::vector<sgNodeID_t> out;
    int degree;
};*/


class CompressionAnalyzer {
public:
    CompressionAnalyzer(SequenceGraph &, uint64_t,  std::string);

    std::vector<NodeCompressions> compressions;
    std::vector<Repeat> FindGraphRepeats();
    void InitializeLib(std::string , std::string , std::string save_to="" );
    void CalculateCompressions();
    void InitializeLibFromDump(std::string );
    void CalculateRepeatCompressions();
    void InitializeCoreGenomeFinder();
    void FindCoreGenome();
private:
    void InitializeKCI();

        void Calculate(NodeCompressions & , std::string mode="analytic");
        SequenceGraph &sg;

        uint64_t max_mem_gb;
    std::string outfile_prefix;
    std::string outfile_csv_name;

    CoreGenomeFinder cgf;
    KmerCompressionIndex kci;
    RepeatAnalyzer ra;
};


#endif //BSG_COMPRESSIONANALYZER_H
