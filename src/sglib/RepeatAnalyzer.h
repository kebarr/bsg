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

int add(int a, int b){
    return a+b;
}

int do_operation(int (*op)(int a, int b) , int arg1, int arg2){
    return op(arg1,arg2);
}

struct RepeatCompressions {

    std::string compression_function_name;
//    double (*compression_function)(sgNodeID_t, uint16_t, int) ;
    double (*compression_function)(sgNodeID_t) ;

    Repeat& repeat;
    sgNodeID_t repeated_contig;

    std::vector<double > in_compressions;
    std::vector<double > out_compressions;
    double repeat_compression;
    //compression function is function used to calculate representation of read set in contig
   // RepeatCompressions(std::string cfn, double (*compression_function)(sgNodeID_t, uint16_t=10, int=0), Repeat& repeat) :
    RepeatCompressions(std::string cfn, double (*compression_function)(sgNodeID_t), Repeat& repeat) :
            compression_function_name(cfn), compression_function(compression_function), repeat(repeat){
        repeat_compression = compression_function(repeat.repeated_contig);
       this->repeated_contig = repeat.repeated_contig;
       for (auto in:repeat.in_contigs){ this->in_compressions.push_back(compression_function(in));}
       ;

       for (auto out:repeat.out_contigs) this->out_compressions.push_back(compression_function(out));
    }

};

class RepeatAnalyzer{
public:
    RepeatAnalyzer(SequenceGraph& _sg, std::string lib_name="");
    void FindRepeats(std::string name_base="rep", int limit=-1, int rep_min = 3000, int in_min=1000,int out_min=1000);

    void OutputRepeats(std::string fname, std::vector<size_t > to_include={});

    Repeat RepeatReduction(Repeat );
    std::vector<Repeat> repeats;

private:
    SequenceGraph & sg;

};

#endif //BSG_REPEATANALYZER_H
