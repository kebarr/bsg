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


struct NodeCompressions {
    std::string lib_name_r1;
    std::string lib_name_r2;
    int index;

    //std::map<sgNodeID_t, double> compressions; - just use sg indexing
    std::vector< double> compressions;
    std::vector<std::vector<sgNodeID_t >> canonical_repeats;
};


class CompressionAnalyzer {
public:
    CompressionAnalyzer(SequenceGraph &, uint64_t,  std::string);

    std::vector<NodeCompressions> compressions;
    void InitializeLib(std::string , std::string , std::string save_to="" );
    void CalculateCompressions();
    void InitializeLibFromDump(std::string );

private:
    void InitializeKCI();

    std::vector<std::vector<double>> AnalyseRepeat(std::vector<std::vector<double>>  repeat_compressions, double tolerance=0.8, double diff_threshold=0.8);
        void Calculate(NodeCompressions & );
    std::vector<double > CompressionStats(std::vector<double> res);
        SequenceGraph &sg;
    uint64_t max_mem_gb;
    std::string outfile_name;
    std::string outfile_csv_name;
    std::string outfile_csv_name1;

    std::string outfile_csv_name2;
    std::string outfile_csv_name3;
    std::string outfile_csv_name4;
    std::string outfile_csv_name5;
    std::string outfile_csv_name6;
    std::string outfile_csv_name7;


    KmerCompressionIndex kci;
};


#endif //BSG_COMPRESSIONANALYZER_H
