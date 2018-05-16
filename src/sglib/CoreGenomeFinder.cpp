//
// Created by Katie Barr (EI) on 15/05/2018.
//

#include "CoreGenomeFinder.h"

CoreGenomeFinder::CoreGenomeFinder(SequenceGraph& _sg, KmerCompressionIndex& _kci, CoreGenomeParams _gcp):sg(_sg), kci(_kci), gcp(_gcp){
};

void CoreGenomeFinder::InitialiseNodeMetrics(){
    std::cout << "Initialising Node Metrics for " << sg.nodes.size() << " nodes"<< std::endl;
    int count =0;
    int candidate_count = 0;
    for (auto node:sg.nodes) {
        std::vector<uint64_t> nkmers;
        StringKMerFactory skf(node.sequence, 31);
        nkmers.reserve(node.sequence.size());
        skf.create_kmers(nkmers);

        std::vector<sgNodeID_t> in_contigs;
        std::vector<sgNodeID_t> out_contigs;
        auto in = sg.get_bw_links(count);
        for (auto c:in) {

            in_contigs.push_back(c.dest);
        }


        auto out = sg.get_fw_links(count);
        for (auto c:out) {

            out_contigs.push_back(c.dest);
        }
        NodeMetrics nm(kci.graph_kmers, nkmers, node, count, gcp, in_contigs, out_contigs);
        nms.push_back(nm);
if (nm.candidate_core){
    candidate_count += 1;
    std::cout << " len : " << nm.sequence_length << " id " << nm.id << " in " << nm.in_contigs.size() << " out " << nm.out_contigs.size() <<
        " kmers " << nm.kmers.size() << " unique " <<  nm.unique_kmer_mask.size() <<  std::endl;
}
        count +=1;
    }
    std::cout << "found " << candidate_count << " candidate core genome nodes " << std::endl;
};

void CoreGenomeFinder::CalculateMetricForReadSet(std::string function_name, double (*compression_function)(sgNodeID_t, KmerCompressionIndex&), int read_set){
    kci.current_lib = read_set;
    for (auto n:nms) {
        if (n.candidate_core) {
            if (n.lib_vals.find(function_name) == n.lib_vals.end()) {
                n.lib_vals[function_name] = {};
            }
            if (n.mapped_libs.find(function_name) == n.mapped_libs.end()) {
                n.mapped_libs[function_name] = {};

            }
            auto res = compression_function(n.id, kci);
            n.lib_vals[function_name].push_back(res);
            if (res > gcp.lib_kmer_thresh) {
                n.mapped_libs[function_name].push_back(read_set);
            }
        }
    }
};

/* this will probably savetine if i can work out how....
void CoreGenomeFinder::DumpNodeMetrics(std::string filename) {
    std::ofstream of(filename);
    //read-to-tag
    uint64_t kcount=kci.graph_kmers.size();
    of.write((const char *) &kcount,sizeof(kcount));
    of.write((const char *) kci.graph_kmers.data(),sizeof(KmerCount)*kcount);
    for (auto m:nms){
        of.write((const char *) m.kmers.data(), sizeof(uint64_t)* m.kmers.size());
        of.write((const char *) m.unique_kmer_mask, sizeof(uint64_t)* m.kmers.size());

    }
    //read-to-node
    uint64_t ccount=read_counts.size();
    of.write((const char *) &ccount,sizeof(ccount));
    if (lib == -1) {
        for (auto i = 0; i < ccount; ++i) {
            of.write((const char *) read_counts[i].data(), sizeof(uint16_t) * kcount);
        }
    } else {
        of.write((const char *) read_counts[lib].data(), sizeof(uint16_t) * kcount);
    }
};*/

