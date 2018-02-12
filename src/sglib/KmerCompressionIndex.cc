//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"



KmerCompressionIndex::KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem): sg(_sg) {
    max_mem = _max_mem;
}

void KmerCompressionIndex::index_graph(){
    const int k = 31;
    const int max_coverage = 1;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, max_mem, 0, max_coverage,
                                                                          output_prefix);



    std::cout << "Indexing graph... " << std::endl;
    graph_kmers = kmerCount_SMR.process_from_memory();

    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
    // [0]- total records generated. [2] - number of reader records, graph treated like fasta

    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << " gk size: " << graph_kmers.size() <<  std::endl;

}

void KmerCompressionIndex::reindex_graph(){
    const int k = 31;
    const int max_coverage = 1;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, max_mem, 0, max_coverage,
                                                                              output_prefix);



    std::cout << "Re-indexing graph... " << std::endl;
    auto new_graph_kmers = kmerCount_SMR.process_from_memory();
    uint64_t deleted=0,changed=0,equal=0;
    //std::sort(new_graph_kmers.begin(),new_graph_kmers.end());
    for (auto i=0,j=0;i<graph_kmers.size() and j<new_graph_kmers.size();++j){
        while (i<graph_kmers.size() and graph_kmers[i].kmer<new_graph_kmers[j].kmer) {
            graph_kmers[i].count=0;
            ++deleted;
            ++i;
        }
        if (i<graph_kmers.size() and graph_kmers[i].kmer==new_graph_kmers[j].kmer){
            if (graph_kmers[i].count==new_graph_kmers[j].count) ++equal;
            else {
                graph_kmers[i].count=new_graph_kmers[j].count;
                ++changed;
            }
            ++i;
        }
    }
    std::cout << deleted << " deleted,   "<<changed<<" changed,   "<<equal<<" equal"<<std::endl;
    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::load_from_disk(std::string filename) {
    std::ifstream inf(filename);
    //read-to-tag
    uint64_t kcount;
    inf.read(( char *) &kcount,sizeof(kcount));
    graph_kmers.resize(kcount);
    inf.read(( char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    //read-to-node
    uint64_t ccount;
    inf.read(( char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        read_counts.emplace_back();
        read_counts.back().resize(kcount);
        inf.read(( char *) read_counts.back().data(), sizeof(uint16_t) * kcount);
    }

}

void KmerCompressionIndex::save_to_disk(std::string filename) {
    std::ofstream of(filename);
    //read-to-tag
    uint64_t kcount=graph_kmers.size();
    of.write((const char *) &kcount,sizeof(kcount));
    of.write((const char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    //read-to-node
    uint64_t ccount=read_counts.size();
    of.write((const char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        of.write((const char *) read_counts[i].data(), sizeof(uint16_t) * kcount);
    }
}

void KmerCompressionIndex::start_new_count(){
    read_counts.emplace_back();
    read_counts.back().resize(graph_kmers.size(),0);
}

void KmerCompressionIndex::add_counts_from_file(std::string filename) {

    std::cout<<"counting from file"<<std::endl;

    FastqReader<FastqRecord> fastqReader({0},filename);
    std::atomic<uint64_t> present(0), absent(0), rp(0);
    std::unordered_map<uint64_t,uint64_t> kmer_map;
    for (uint64_t i=0;i<graph_kmers.size();++i) kmer_map[graph_kmers[i].kmer]=i;

#pragma omp parallel shared(fastqReader)
    {
        std::vector<uint16_t> thread_counts(read_counts.back().size(),0);
        FastqRecord read;
        std::vector<KmerCount> readkmers;
        KmerCountFactory<FastqRecord> kf({31});

        bool c;
#pragma omp critical(fastqreader)
        c = fastqReader.next_record(read);
        while (c) {
            readkmers.clear();
            //process tag if 10x! this way even ummaped reads get tags
            kf.setFileRecord(read);
            kf.next_element(readkmers);

            for (auto &rk:readkmers) {
                /*auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), rk);
                if (nk != graph_kmers.end() and nk->kmer == rk.kmer) {
                    auto offset = nk - graph_kmers.begin();
                    ++thread_counts[offset];
                    ++present;
                } else ++absent;*/
                auto findk = kmer_map.find(rk.kmer);
                if (kmer_map.end() != findk){
                    ++thread_counts[findk->second];
                    ++present;
                } else ++absent;


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
            for (uint64_t i=0;i<read_counts.back().size();++i) read_counts.back()[i]+=thread_counts[i];
        }
    }
    std::cout << rp << " reads processed "<< present <<" / " << present+absent << " kmers found" << std::endl;
}



void KmerCompressionIndex::compute_compression_stats(size_t lib) {
    //compute mean, median and mode, as of now, only use the first read count

    uint64_t covuniq[1001];
    for (auto &c:covuniq)c=0;
    uint64_t tuniq=0,cuniq=0;
    for (uint64_t i=0; i<graph_kmers.size(); ++i){
        if (graph_kmers[i].count==1){
            tuniq+=read_counts[lib][i];
            ++cuniq;
            ++covuniq[(read_counts[lib][i]<1000 ? read_counts[lib][i] : 1000 )];
        }
    }
    uint64_t cseen=0,median=0;
    while (cseen<cuniq/2) {cseen+=covuniq[median];++median;};
    uint64_t mode=0;
    for (auto i=0;i<1000;++i) if (covuniq[i]>covuniq[mode]) mode=i;
    std::cout << "Mean coverage for unique kmers:   " << ((double)tuniq)/cuniq <<std::endl;
    std::cout << "Median coverage for unique kmers: " << median <<std::endl;
    std::cout << "Mode coverage for unique kmers:   " << mode <<std::endl;

    if (median<.9*mode or median>.9*mode ) std::cout<<"WARNING -> median and mode highly divergent"<<std::endl;
    uniq_mode=mode != 0 ? mode:1;// this was usually 0... makes no sense to have a kmer that appears 0 times!!

}

void KmerCompressionIndex::dump_histogram(std::string filename, uint16_t dataset) {
    std::ofstream kchf(filename);
    uint64_t covuniq[1001];
    for (auto &c:covuniq)c=0;
    uint64_t tuniq=0,cuniq=0;
    for (uint64_t i=0; i<graph_kmers.size(); ++i){
        if (graph_kmers[i].count==1){
            tuniq+=read_counts[dataset][i];
            ++cuniq;
            ++covuniq[(read_counts[dataset][i]<1000 ? read_counts[dataset][i] : 1000 )];
        }
    }
    for (auto i=0;i<1000;++i) kchf<<i<<","<<covuniq[i]<<std::endl;
}

double KmerCompressionIndex::compute_compression_for_node(sgNodeID_t _node, uint16_t max_graph_freq, uint16_t dataset) {

    auto & node=sg.nodes[_node>0 ? _node:-_node];

    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence,31);
    skf.create_kmers(nkmers);
    //std::cout << node.sequence << std::endl;
    int counter = 0;
    int counter_pres = 0;

    uint64_t kcount=0,kcov=0;
    std::cout << "kmers in nodes: " << nkmers.size() << std::endl;
    for (auto &kmer : nkmers){
        // find kmer in graph kmer with count > 0?
        auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (nk!=graph_kmers.end() and nk->kmer == kmer and nk->count==1) {
            counter +=1;
            ++kcount;// inrement number of (unique?? ) kmers on node
            kcov+=read_counts[dataset][nk-graph_kmers.begin()]; // inrement coverage by count for this kmer in read set
        } else if (nk!=graph_kmers.end() and nk->kmer == kmer and nk->count> 1) {
            counter_pres+= 1;
        }

    }
    std::cout << counter << " kmers found " << counter_pres << " kmers  count > 1 \n";
    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return (((double) kcov)/kcount )/uniq_mode;
}