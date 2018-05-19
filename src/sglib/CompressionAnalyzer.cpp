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
    kci.add_counts_from_file(lib_name_r1);
    kci.add_counts_from_file(lib_name_r2);
    kci.compute_compression_stats(kci.read_counts.size()-1);

    if (save_to != ""){
        kci.save_to_disk(save_to, nc.index);
    }

};

void CompressionAnalyzer::InitializeCoreGenomeFinder( ){
    cgf.InitialiseNodeMetrics();

};


double compute_kmers_present_for_node(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci, int dataset){

    int counter = 0;

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        counter +=1;
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count > 0 && kci.read_counts[dataset][kci.kmer_map[kmer]] > 0 ){// should scale non uniwue kmers by number occurnces in graph
            kcov+=1;}

    }
    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov/(double)nkmers.size();
};


double evaluate_contiguous_matches(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci, int dataset){
    //of the matching kmers- how many map contiguously
};

double unique_kmers_present(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci, int dataset) {

    int counter = 0;

    uint64_t kcount = 0, kcov = 0;
    for (auto &kmer : nkmers) {
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count ==1) {
            counter += 1;
            if (kci.read_counts[dataset][kci.kmer_map[kmer]] >
                                         0) {// should scale non uniwue kmers by number occurnces in graph
                kcov += 1;
            }

        }

    }       // number of times unique kmers in this node appear in reads, scaled by total unique kmers
        return kcov / (double) counter;
};


double compute_unique_kmers_for_node(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci, int dataset){

    int counter = 0;

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count ==1) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            kcov+=kci.read_counts[dataset][kci.kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }

    // number of times unique kmers in this node appear in reads, scaled by total unique kmers
    return kcov/((double)counter*kci.mean);
};

double read_coverage(std::vector<uint64_t> nkmers, KmerCompressionIndex& kci, int datase){

    int counter = 0;

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (kci.graph_kmers[kci.kmer_map[kmer]].count > 0) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            kcov+=kci.read_counts[datase][kci.kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }
    // count con
    //tig as repesented in read set if more than 70%? kmers present
    // this will end up same as one of others...
    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov/((double)nkmers.size())*kci.mean;
};


void CompressionAnalyzer::FindCoreGenome() {

CoreGenomeParams gcp;
    gcp = CoreGenomeParams();
    std::cout << " finding core based on " << kci.read_counts.size() << " varieties " << std::endl;
    std::map<std::string , double (*)(std::vector<uint64_t> , KmerCompressionIndex&, int)> metric_map ={{"unique_kmers_present",unique_kmers_present}, {"read_coverage", read_coverage},{"compute_kmers_present_for_node", compute_kmers_present_for_node}, {"compute_unique_kmers_for_node", compute_unique_kmers_for_node}};
    InitializeCoreGenomeFinder();
    for (auto m:metric_map) {
 auto res =       cgf.EvaluateMetric(m.first, m.second);

        std::ofstream outfile;
        outfile.open(outfile_prefix + m.first + ".csv");
        for (auto n:cgf.nms){
            outfile << sg.oldnames[n.id] << ", " ;
        }
        outfile << std::endl;
        for (auto r:res){
            for (auto v:r){
                outfile << v << ", ";
            }
            outfile << std::endl;

        }
    }
    cgf.SelectCoreGenome();
    cgf.OutputNodeMetrics(outfile_prefix + "metrics.txt");
    cgf.OutputCoreFasta(outfile_prefix );

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

};

void  CompressionAnalyzer::UniqueContigsForLibs(std::vector<int> libs, int lib){




// to find contigsunique to each read set, need kmers unique to each, contigs with > some % unique kmers are considered unique to that variety

void  CompressionAnalyzer::UniqueKmersForLib(std::vector<int> libs, int lib){

};

// as need unique contigs... just getgraph kmers
std::set<uint64_t > CompressionAnalyzer::count_kmers_from_file(std::string filename, int max_freq) {


    FastqReader<FastqRecord> fastqReader({0},filename);
    std::atomic<uint64_t> present(0), absent(0), rp(0);
    std::vector<uint16_t > read_count;
    size_t  sum = 0;
    size_t  kmers_in_reads = 0;
    std::unordered_map<uint64_t,uint64_t> kmer_map;

#pragma omp parallel shared(fastqReader)
    {
        std::vector<uint16_t> thread_counts(read_count,0);
        FastqRecord read;
        std::vector<KmerCount> readkmers;
        KmerCountFactory<FastqRecord> kf({31});

        bool c;
#pragma omp critical(fastqreader)
        c = fastqReader.next_record(read);
        while (c) {
            readkmers.clear();
            kf.setFileRecord(read);
            kf.next_element(readkmers);

            for (auto &rk:readkmers) {

                auto findk = kmer_map.find(rk.kmer);
                if (kmer_map.end() != findk){
                    ++thread_counts[findk->second];
                } else {
                    kmer_map[rk.kmer] =(int)thread_counts.size();
                    thread_counts.push_back(1);
                }


            }
            uint64_t a=++rp;
            if (a % 100000 == 0)
                std::cout << rp << " reads processed " << present << " / " << present + absent << " kmers found"
                          << std::endl;
#pragma omp critical(fastqreader)
            c = fastqReader.next_record(read);
        }

        //Update shared counts
#pragma omp critical(countupdate)
        {

            for (uint64_t i=0;i<read_count.size();++i) {
                read_count[i]+=thread_counts[i];
                sum += thread_counts[i];
                if (thread_counts[i] != 0){
                    kmers_in_reads += 1;
                }
            }
        }
    }
    std::set<uint64_t > kmers;
    for (auto k:kmer_map){
        auto count = read_count[k.second];
        if (count < max_freq){
            kmers.insert(k.first);
        }
    }

    std::cout << rp << " reads processemers, kmers  in reads: " << kmers.size() << std::endl;

    return; kmers.
}
