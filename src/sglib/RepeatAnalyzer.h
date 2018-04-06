//
// Created by Katie Barr (EI) on 06/04/2018.
//

#ifndef BSG_REPEATANALYZER_H
#define BSG_REPEATANALYZER_H

#include <string>
#include <map>
#include "sglib/SequenceGraph.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <sglib/KmerCompressionIndex.hpp>

struct Repeat {
    std::string name_base;
    sgNodeID_t repeated_contig; // theactual node can always be retrieved from sg
    std::vector<sgNodeID_t> in_contigs;
    std::vector<sgNodeID_t> out_contigs;

};

class RepeatAnalyzer{
    RepeatAnalyzer(SequenceGraph &_sg, std::string lib_name="");
    void FindRepeats(std::string name_base="rep", int limit=-1, int rep_min = 3000, int in_min=1000,int out_min=1000);

    void OutputRepeats(std::string fname, std::vector<size_t > to_include={});

    void RepeatReduction(std::vector<sgNodeID_t > nodes );

private:
    SequenceGraph & sg;
std::vector<Repeat> repeats;
};

#endif //BSG_REPEATANALYZER_H
