//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"
#include "RepeatAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_name(outfile_prefix + ".txt"), outfile_csv_name(outfile_prefix+".csv"),kci(this->sg, this->max_mem_gb*1024L*1024L*1024L),
                                                                                                               ra(this->sg, this->kci), cgf(this->sg, this->kci, CoreGenomeParams()){

    InitializeKCI();

};

void CompressionAnalyzer::InitializeLibFromDump(std::string lib_name) {
    std::cout << "Retrieving kci of " << lib_name << std::endl;
    NodeCompressions nc;
    nc.lib_name_r1 = lib_name;
    nc.lib_name_r2 = lib_name;
    nc.index= static_cast<int>(compressions.size());
    std::cout << "loaded kci dump from file, " << lib_name << std::endl;
    kci.load_from_disk(lib_name);
    std::cout << "loaded kci from file, " << kci.read_counts.size() << " read set dumps loaded" << std::endl;

}

void CompressionAnalyzer::InitializeLib(std::string lib_name_r1, std::string lib_name_r2, std::string save_to="") {
    std::cout << "Initializing  compression analysis of " << lib_name_r1 << " to save to: " << save_to << std::endl;

    NodeCompressions nc = {lib_name_r1, lib_name_r2, static_cast<int>(compressions.size())};
    //nc.index = compressions.size();
    nc.compressions.resize(sg.nodes.size(), -1);
    compressions.push_back(nc);
    kci.start_new_count();

    std::cout << sg.nodes.size()<< std::endl;
    std::cout <<kci.sg.nodes.size()<< std::endl;
    kci.add_counts_from_file(lib_name_r1);
    kci.add_counts_from_file(lib_name_r2);
    if (save_to != ""){
        kci.save_to_disk(save_to, nc.index);
    }

};

void CompressionAnalyzer::InitializeCoreGenomeFinder(){
    cgf.InitialiseNodeMetrics();

};


double compute_kcov_for_node(sgNodeID_t _node, KmerCompressionIndex& kci){
    auto & node=kci.sg.nodes[_node>0 ? _node:-_node];
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence,31);
    skf.create_kmers(nkmers);
    //std::cout << node.sequence << std::endl;
    int counter = 0;

    uint64_t kcount=0,kcov=0;
    //std::cout << "kmers in node: "<< _node << ", " << nkmers.size() << std::endl;
    for (auto &kmer : nkmers){
        // find kmer in graph kmer with count > 0?
        // need index of kmer in hraph_kmera - must be a better eay

        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count > 0) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            ++kcount;// inrement number of (unique??- now removed count = 1 ) kmers on node
            kcov+=kci.read_counts[kci.current_lib][kci.kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }
    //std::cout << "kcount: " << kcount << " kcov " << kcov << " kcountcount: " << kcountcount << " kcountu: " << kcountu << " kcovu " << kcovu << " kcountcountu: " << kcountcountu <<std::endl;

    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov;
};

void CompressionAnalyzer::FindCoreGenome() {
    std::cout << " finding core based on " << kci.read_counts.size() << " varieties " << std::endl;
    InitializeCoreGenomeFinder();

        for (int i = 0; i < kci.read_counts.size(); i++) {
            cgf.CalculateMetricForReadSet("basic", compute_kcov_for_node,  i);

        }

};


void CompressionAnalyzer::InitializeKCI(){
    std::cout << "Initializing kmer copression index, indexig sequene graph" << std::endl;
    //KmerCompressionIndex kci(sg, max_mem_gb * 1024L * 1024L * 1024L);
    //this->kci = kci;
    kci.index_graph();

    std::cout << sg.nodes.size()<< std::endl;
    std::cout <<kci.sg.nodes.size()<< std::endl;
    std::cout << "Initialization complete" << std::endl;

     std::ofstream outfile_csv;

    outfile_csv.open(outfile_csv_name, std::ofstream::out );

    for (size_t counter = 0; counter < sg.nodes.size(); counter++) {
        outfile_csv << sg.oldnames[counter] << ", ";


    }
    outfile_csv << std::endl;

    if (kci.read_counts.size()>0) {
        kci.compute_compression_stats();
        kci.dump_histogram(outfile_name + "kci_histogram.csv");
    }

}


void CompressionAnalyzer::CalculateRepeatCompressions() {
   int resolved_repeat_indices = 0;
for (int i=0;i <kci.read_counts.size(); i++) {
    kci.current_lib = i;
   auto c= ra.compressions_for_read_set(compute_kcov_for_node);
    this->compressions[i].compressions = c;
}
// now have all compressions for each lib (though not sure its needed), and each repeat has the compressions from all libs so that we can see whether the different varieties resolve it

};
