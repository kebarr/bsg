//
// Created by Katie Barr (EI) on 15/05/2018.
//

#ifndef BSG_COREGENOMEFINDER_H
#define BSG_COREGENOMEFINDER_H


#include <string>
#include <map>
#include "sglib/SequenceGraph.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <sglib/KmerCompressionIndex.hpp>



struct CoreGenomeParams {
    int sequence_length_thresh=1000;
    double kmer_thresh=0.1; // proportion ofunique  kmers requireed in contigs
    int min_links=1; // miimum number of contigs are core genome candidate shpuld map to
    int min_libs=5;// min varieties that candidates should be present in
    double lib_kmer_thresh = 0.8; // how many kmers in cntig shouls be present in read set to consider that contig included by that read set
};

struct NodeMetrics{
    Node & node;
    CoreGenomeParams gcp;
    //sav typing .node all the time
    size_t  sequence_length;
    sgNodeID_t  id;
    std::vector<uint64_t > kmers;
    std::vector<bool> unique_kmer_mask ;
    double kmer_thresh;
    std::vector<sgNodeID_t> in_contigs;
    std::vector<sgNodeID_t> out_contigs;
    bool candidate_core;
    bool sane_flow;
    std::map<std::string,  std::vector<int>> mapped_libs; // for each metric, libs that met threshold
    //string is name of metrics function and vector is result of that function each read set
    std::map<std::string, std::vector<double> > lib_vals;
    NodeMetrics(std::vector<KmerCount> & graph_kmers, std::vector<uint64_t> nkmers, Node & _node, sgNodeID_t _id, CoreGenomeParams _gcp, std::vector<sgNodeID_t> _in_contigs, std::vector<sgNodeID_t> _out_contigs):
            node(_node), id(_id), sequence_length(this->node.sequence.size()), gcp(_gcp), in_contigs(_in_contigs), out_contigs(_out_contigs){
        unique_kmer_mask.resize(nkmers.size());
        int counter = 0;
        int i = 0;
        for (auto &kmer : nkmers){
            // find kmer in graph kmer with count > 0?
            auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
            if (nk!=graph_kmers.end() and nk->kmer == kmer and nk-> count == 1) {
                counter +=1;
                unique_kmer_mask[i] = true;

            }

        }
        double unique_kmers = counter/(double)nkmers.size();
        if (sequence_length > gcp.sequence_length_thresh && in_contigs.size() + out_contigs.size() > gcp.min_links && unique_kmers > gcp.kmer_thresh ){
            candidate_core = true;
        }
    }


};

class CoreGenomeFinder {
public:
    CoreGenomeFinder(SequenceGraph& _sg, KmerCompressionIndex& _kci, CoreGenomeParams _gcp);
    std::vector<NodeMetrics > nms;
    void DumpNodeMetrics(std::string );
    void InitialiseNodeMetrics();
    void CalculateMetricForReadSet(std::string function_name, double (*compression_function)(sgNodeID_t, KmerCompressionIndex&), int read_set);

private:

    SequenceGraph & sg;
    KmerCompressionIndex & kci;
    CoreGenomeParams gcp;
};


#endif //BSG_COREGENOMEFINDER_H
