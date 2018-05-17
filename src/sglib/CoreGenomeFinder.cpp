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
    int seq_len = 0;
    int candidate_seq_len = 0;
    for (auto node:sg.nodes) {
        std::vector<uint64_t> nkmers;
        StringKMerFactory skf(node.sequence, 31);
        seq_len += node.sequence.size();
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
        NodeMetrics nm(kci, nkmers, node, count, gcp, in_contigs, out_contigs);
        nms.push_back( nm);
if (nm.candidate_core){
    candidates.push_back(nm.id);
    candidate_count += 1;
    candidate_seq_len += node.sequence.size();
    std::cout << " len : " << nm.sequence_length << " id " << nm.id << " in " << nm.in_contigs.size() << " out " << nm.out_contigs.size() <<
        " kmers " << nm.kmers.size() <<  std::endl;
}
        count +=1;

    }
    this->candidates = candidates;
    std::cout << "found " << candidate_count << " candidate core genome nodes with  seq len " << candidate_seq_len << " out of " << seq_len << " bases in graph" << std::endl;
};

std::vector<double > CoreGenomeFinder::CalculateMetricForReadSet(std::string function_name, double (*compression_function)(std::vector<uint64_t> , KmerCompressionIndex&), std::string filename,int read_set){
    kci.current_lib = read_set;
    std::ofstream outfile;
    outfile.open(filename);
    int mapped = 0;
    std::vector<double > vals;
    std::cout << " calculating " << function_name << " for " << read_set << std::endl;
    for (auto i:candidates) {
        auto n = nms[i];
        outfile << n.id << ", ";

            auto res = compression_function(n.kmers, kci);
        vals.push_back(res);
        nms[i].update_metric(function_name,res);
            if (res > gcp.lib_kmer_thresh) {

                nms[i].increment_number_libs_mapped(function_name);
                mapped += 1;
                std::cout <<  " libs mapped: " << nms[i].number_libs_mapped[function_name] << std::endl;

            }
            std::cout << "node " << n.id << "res " << res << " libs mapped: " << nms[i].number_libs_mapped[function_name] << std::endl;
        //this->nms[i] = n;

    }
outfile<< std::endl;
    for (auto r:vals){
        outfile << r << ", " ;
    }
    outfile<< std::endl;

    std::cout << " read_set " << read_set << " mapped to " << mapped << " candidate contigs" << std::endl;
    return  vals;
};


void CoreGenomeFinder::OutputNodeMetrics(std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto n:nms){
        outfile << "id: " <<  n.id << " length, " << n.sequence_length<< " core candidate " << n.candidate_core << " core: " << n.core << " ";
        for (auto m:n.lib_vals) {
            outfile << m.first << ": ";
            for (auto v: m.second) {
                outfile << v << ", ";
            }
            outfile << std::endl;
        }
    }
}

void CoreGenomeFinder::OutputCoreFasta(std::string filename){
    std::ofstream outfile;
    outfile.open(filename);
    std::cout <<  " output contigs to " << filename << std::endl;

    for (auto n:core){
        outfile << ">" << sg.oldnames[nms[n].id] << std::endl << sg.nodes[n].sequence  << std::endl;

    }
};

void CoreGenomeFinder::SelectCoreGenome(){
    std::cout << " choosing core contigs " << std::endl;
    int len=0;
    for (auto i:candidates) {
        auto n = nms[i];
        for (auto t: n.number_libs_mapped) {
            if (t.second > gcp.min_libs) {
                core.push_back(n.id);
                len += n.sequence_length;
                this->nms[i].core = true;
                break;
            }
        }
    }
    std::cout << core.size() << " contigs cound with total length " << len << std::endl;

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

