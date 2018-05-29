//
// Created by Katie Barr (EI) on 15/05/2018.
//

#include "CoreGenomeFinder.h"

CoreGenomeFinder::CoreGenomeFinder(SequenceGraph& _sg, KmerCompressionIndex& _kci, CoreGenomeParams _gcp):sg(_sg), kci(_kci), gcp(_gcp){
};

void CoreGenomeFinder::InitialiseNodeMetrics() {
    std::cout << "Initialising Node Metrics for " << sg.nodes.size() << " nodes" << std::endl;
    int count = 0;
   // int seq_len = 0;

#pragma omp parallel for
    for (int i = 0 ; i <sg.nodes.size() ; i++) {
        auto node = sg.nodes[i];
        if (count % 10000 == 0){
            std::cout << "Initialized metrics for " << std::to_string(count) << " nodes, of which " <<  std::to_string(count) << " are candidates for core genomes"<< std::endl;
        }
        std::vector<uint64_t> nkmers;
        StringKMerFactory skf(node.sequence, 31);
        //seq_len += node.sequence.size();
        nkmers.reserve(node.sequence.size());
        skf.create_kmers(nkmers);

        std::vector<sgNodeID_t> in_contigs;
        std::vector<sgNodeID_t> out_contigs;
        /*auto in = sg.get_bw_links(count);
        for (auto c:in) {

            in_contigs.push_back(c.dest);
        }


        auto out = sg.get_fw_links(count);
        for (auto c:out) {

            out_contigs.push_back(c.dest);
        }*/
        NodeMetrics nm(kci, nkmers, node, count, gcp, in_contigs, out_contigs);
#pragma omp critical(nm)
{
        nms.push_back(nm);
            candidates.push_back(nm.id);

        count += 1;
    }

    };
        this->candidates = candidates;
        std::cout << "found " << count << " candidate core genome nodes with  seq len " <<
                   std::endl;

};



std::vector<std::vector<double > > CoreGenomeFinder::EvaluateMetric( std::string  function_name, double (*compression_function)(std::vector<uint64_t> , KmerCompressionIndex&, int)){
    auto metrics_so_far = function_names.size();
    std::vector<std::vector<double> > res;
    function_names[function_name] = metrics_so_far;
    std::cout << "evaluating  " << function_name << metrics_so_far << "th metric "  << std::endl;
    for (int i = 0; i < kci.read_counts.size(); i++) {
        auto r = CalculateMetricForReadSet(function_name, compression_function, i);
        res.push_back(r);

    }


    return  res;
};


std::vector<double > CoreGenomeFinder::CalculateMetricForReadSet(std::string function_name, double (*compression_function)(std::vector<uint64_t> , KmerCompressionIndex&, int),
                                                                 int read_set){

    int mapped = 0;
    std::vector<double > vals;
    std::cout << " calculating " << function_name << " for " << read_set << std::endl;
    for (auto i:candidates) {
        auto n = nms[i];
        nms[i].add_metric(function_names[function_name]);
            auto res = compression_function(n.kmers, kci, read_set);
        vals.push_back(res);

        nms[i].update_metric(function_names[function_name],res);
            if (res > gcp.lib_thresh[function_names[function_name]]) {

                nms[i].increment_number_libs_mapped(function_names[function_name]);
                mapped += 1;

            }
        //this->nms[i] = n;

    }

    std::cout << " read_set " << read_set << " mapped to " << mapped << " candidate contigs with " << function_name  << std::endl;
    return  vals;
};


void CoreGenomeFinder::OutputNodeMetrics(std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto n:nms){
        outfile << " id: " <<  sg.oldnames[n.id] << " length, " << n.sequence_length<< " core candidate " << n.candidate_core  << " ";
        outfile << std::endl;
        for (int i=0; i < n.lib_vals.size() ; i++) {
            outfile << " metric " << i << ": ";
            for (auto v: n.lib_vals[i]){
                outfile << v << ", ";
            }
            outfile << std::endl;
        }
    }
}

void  CoreGenomeFinder::OutputCoreFasta(std::string filename){
    for (auto m: core){
        OutputCoreFastaForMetric(filename, m.first);
    }
}

void CoreGenomeFinder::OutputCoreFastaForMetric(std::string filename, std::string const  metric){
    std::ofstream outfile;
    outfile.open(filename+metric+".fasta");
    std::cout <<  " output contigs to " << filename+metric+".fasta" << std::endl;

    for (auto n:core[metric]){
        outfile << ">" << sg.oldnames[nms[n].id] << std::endl << sg.nodes[n].sequence  << std::endl;

    }
};

void CoreGenomeFinder::SelectCoreGenome(){
    std::cout << " choosing core contigs " << std::endl;
    std::vector<int> len;
    int count = 0;
    for (auto fn:function_names){
        len.push_back(0);

        for (auto i:candidates) {
            auto n = nms[i];


            // contig is considered core according to some metric if enough other varieties contain in
            if (n.number_libs_mapped[fn.second] > gcp.min_libs) {
                core[fn.first].push_back(n.id);
                len[count] += n.sequence_length;
                nms[i].is_core_for_metric(fn.second);
            }
        }
            count += 1;
        std::cout  <<  core[fn.first].size() << " contigs cound with total length " << len[fn.second] << " for " << fn.first << std::endl;

    }};

void CoreGenomeFinder::AlternateParams(std::vector<CoreGenomeParams > cgp_new){

    for (auto p:cgp_new){
        std::vector<NodeMetrics  >  nms_new;
        for (auto n:nms){
            if (n.sequence_length > p.sequence_length_thresh && n.in_contigs.size() + n.out_contigs.size() >= gcp.min_links && std::count(n.unique_kmer_mask.begin(), n.unique_kmer_mask.end(), true)/ (double) n.unique_kmer_mask.size() > p.kmer_thresh ){
                nms_new.push_back(n);
            }

        }
    }
;}

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

