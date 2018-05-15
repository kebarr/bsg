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


//compressions for each lib should be stored together
struct RepeatCompressions {

    std::string compression_function_name;
//    double (*compression_function)(sgNodeID_t, uint16_t, int) ;
    double (*compression_function)(sgNodeID_t, KmerCompressionIndex&) ;

    sgNodeID_t repeated_contig;

    std::vector<double > in_compressions;
    std::vector<double > out_compressions;
    bool in_out_agree= false;
    bool repeat_coverage_sane= false;
    bool repeat_contained= false;
    double repeat_compression;
    //compression function is function used to calculate representation of read set in contig
   // RepeatCompressions(std::string cfn, double (*compression_function)(sgNodeID_t, uint16_t=10, int=0), Repeat& repeat) :
    RepeatCompressions(std::string cfn, double (*compression_function)(sgNodeID_t, KmerCompressionIndex&), sgNodeID_t repeated_contig, KmerCompressionIndex &kci, std::vector<sgNodeID_t > in_contigs, std::vector<sgNodeID_t > out_contigs) :
            compression_function_name(cfn), compression_function(compression_function), repeated_contig(repeated_contig){
        this->repeat_compression = compression_function(repeated_contig, kci);
       double in_c = 0;
       for (auto in:in_contigs){auto res=compression_function(in, kci); this->in_compressions.push_back(res); in_c += res;}
       ;
        double out_c = 0;

       for (auto out:out_contigs){auto res=compression_function(out, kci); this->out_compressions.push_back(res); out_c += res;}
        if(repeat_compression > 0 && in_c >0 && out_c > 0) {
            repeat_contained=true;
            if (std::abs(in_c - out_c) < 0.1 * in_c &&
                std::abs(in_c - out_c) < 0.1 * out_c) { this->in_out_agree = true; }
            if (std::abs(in_c - repeat_compression) < 0.1 * in_c &&
                std::abs(out_c - repeat_compression) < 0.1 * out_c) { this->repeat_coverage_sane = true; }
        }
    }


};


struct Repeat {
    std::string name_base;
    sgNodeID_t repeated_contig; // theactual node can always be retrieved from sg
    std::vector<sgNodeID_t> in_contigs;
    std::vector<sgNodeID_t> out_contigs;
    std::vector<std::pair<int, uint64_t>> reduced_in_contigs;
    std::vector<std::pair<int, uint64_t>> reduced_out_contigs;
    //string is name of compressipon function and vector is result of that function each read set
    std::map<std::string, std::vector<RepeatCompressions> >rc;

};


class RepeatAnalyzer{
public:
    RepeatAnalyzer(SequenceGraph& _sg, KmerCompressionIndex& _kci, std::string lib_name="");
    void FindRepeats(std::string name_base="rep", int limit=-1, int rep_min = 3000, int in_min=1000,int out_min=1000);

    void OutputRepeats(std::string fname, std::vector<size_t > to_include={});

    void RepeatReduction(Repeat );
    std::vector<Repeat> repeats;
    std::vector<double> compressions_for_read_set( double (*compression_function)(sgNodeID_t, KmerCompressionIndex &));

private:
    SequenceGraph & sg;
    KmerCompressionIndex & kci;
    std::vector<std::vector<int > > read_sets_mapping_to_repeats;
    std::vector<std::vector<uint64_t> > ConvertParallelContigsDistinctKmers(std::vector<sgNodeID_t > contigs);


    };

#endif //BSG_REPEATANALYZER_H
