//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"
#include "RepeatAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_prefix(outfile_prefix), outfile_csv_name(outfile_prefix+".csv"),kci(this->sg, this->max_mem_gb*1024L*1024L*1024L),
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


double compute_kcov_for_node(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci){

    int counter = 0;

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count > 0) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            kcov+=kci.read_counts[kci.current_lib][kci.kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }
    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov/(double)nkmers.size();
};


double compute_unique_kmers_for_node(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci){

    int counter = 0;

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count ==1) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            kcov+=kci.read_counts[kci.current_lib][kci.kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }

    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov/(double)nkmers.size();
};


void CompressionAnalyzer::FindCoreGenome() {
    std::vector<std::vector<double > > res1;
    std::vector<std::vector<double > > res2;
    res1.resize(kci.read_counts.size());
    res2.resize(kci.read_counts.size());

    std::cout << " finding core based on " << kci.read_counts.size() << " varieties " << std::endl;
    InitializeCoreGenomeFinder();

        for (int i = 0; i < kci.read_counts.size(); i++) {
            auto r1 = cgf.CalculateMetricForReadSet("compute_kcov_for_node",  compute_kcov_for_node, outfile_prefix + std::to_string(i) + "coverage1.csv", i);
            auto r2 = cgf.CalculateMetricForReadSet("compute_unique_kmers_for_node", compute_unique_kmers_for_node, outfile_prefix + std::to_string(i) + "coverage2`.csv", i);
            res1[i] = r1;
            res2[i] = r2;

        }
    cgf.SelectCoreGenome();
    cgf.OutputNodeMetrics(outfile_prefix + "metrics.txt");
    cgf.OutputCoreFasta(outfile_prefix + "core.fasta");

};


void CompressionAnalyzer::InitializeKCI(){
    std::cout << "Initializing kmer copression index, indexig sequene graph" << std::endl;
    //KmerCompressionIndex kci(sg, max_mem_gb * 1024L * 1024L * 1024L);
    //this->kci = kci;
    kci.index_graph();
    int uniqur = 0;
    for (auto k: kci.graph_kmers){
        if (k.count == 1){
            uniqur += 1;
        }
    }
    std::cout << uniqur << " of " << kci.graph_kmers.size() << " uniwue" << std::endl;

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
        kci.dump_histogram(outfile_prefix + "kci_histogram.csv");
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
