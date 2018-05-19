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
    int min_links=1; // miimum number of vontigs are core genome candidate shpuld connect to
    int min_libs=1;// min varieties that candidates should be present in
    std::vector<double> lib_thresh = {0.8, 0.8}; // thresholds for each metric for a  contig to be considered  included by that read set
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
    bool candidate_core= false;
    std::vector<bool> core; // should just use list with indexes corresponding to fn name
    std::vector<int> number_libs_mapped;
    bool sane_flow;
    //std::map<std::string,  std::vector<int>> mapped_libs; // for each metric, libs that met threshold
    //string is name of metrics function and vector is result of that function each read set
    std::vector<std::vector<double> > lib_vals;
    NodeMetrics(KmerCompressionIndex kci, std::vector<uint64_t> nkmers, Node & _node, sgNodeID_t _id, CoreGenomeParams _gcp, std::vector<sgNodeID_t> _in_contigs, std::vector<sgNodeID_t> _out_contigs):
            node(_node), id(_id), kmers(nkmers), sequence_length(this->node.sequence.size()), gcp(_gcp), in_contigs(_in_contigs), out_contigs(_out_contigs){
        unique_kmer_mask.resize(nkmers.size());
        int counter = 0;
        int i = 0;
        for (auto &kmer : nkmers){
            // find kmer in graph kmer with count > 0?
            //auto nk = std::lower_bound(kci.graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
            if (kci.graph_kmers[kci.kmer_map[kmer]].count == 1){
            //if (nk!=graph_kmers.end() and nk->kmer == kmer and nk-> count == 1){
            //if (graph_kmers[i].count == 1){
            //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
            //if (nk!=graph_kmers.end() and nk->kmer == kmer and nk-> count == 1) {
                counter +=1;
                unique_kmer_mask[i] = true;

            }
            i +=1;


        }
        double unique_kmers = counter/(double)nkmers.size();
        if (sequence_length > gcp.sequence_length_thresh && in_contigs.size() + out_contigs.size() >= gcp.min_links && unique_kmers > gcp.kmer_thresh ){
            candidate_core = true;
        }
    }

    void add_metric(int i){
        if (number_libs_mapped.size() == i) {
            number_libs_mapped.push_back(0);
            lib_vals.emplace_back();
            core.push_back(false);
        }
    }
    void increment_number_libs_mapped(int fn){

        number_libs_mapped[fn] += 1;
    }

    void update_metric(int fn, double value){

            lib_vals[fn].push_back(value);
    }

    void is_core_for_metric(int fn){
        core[fn] = true;
    }
// https://openclassrooms.com/forum/sujet/c-11-use-of-deleted-function no iea what i\m doing!
    //NodeMetrics(NodeMetrics const &);
   // NodeMetrics & operator = (NodeMetrics const &);
};


class CoreGenomeFinder {
public:
    CoreGenomeFinder(SequenceGraph& _sg, KmerCompressionIndex& _kci, CoreGenomeParams _gcp);
    void DumpNodeMetrics(std::string );
    void SelectCoreGenome();
    void InitialiseNodeMetrics();
    void OutputNodeMetrics(std::string  filename);
    void OutputCoreFasta(std::string );
    void EvaluateAllMetrics();

    std::vector<double > CalculateMetricForReadSet(std::string function_name, double (*compression_function)(std::vector<uint64_t> , KmerCompressionIndex&, int),  int );
    std::vector<NodeMetrics  >  nms;
    std::vector< std::vector<double > > EvaluateMetric( std::string , double (*)(std::vector<uint64_t> , KmerCompressionIndex&, int));
    void AlternateParams(std::vector<CoreGenomeParams >);\


        private:
    std::map<std::string, int> function_names;

    void OutputCoreFastaForMetric(std::string , std::string  const );

    std::vector<sgNodeID_t > candidates;
    std::map< std::string, std::vector<sgNodeID_t > >core;

    SequenceGraph & sg;
    KmerCompressionIndex & kci;
    CoreGenomeParams gcp;
    CoreGenomeFinder(CoreGenomeFinder const &);
    CoreGenomeFinder & operator = (CoreGenomeFinder const &);
};


#endif //BSG_COREGENOMEFINDER_H
