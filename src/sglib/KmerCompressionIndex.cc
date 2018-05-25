//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"



KmerCompressionIndex::KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem): sg(_sg) {
    max_mem = _max_mem;
}



void KmerCompressionIndex::index_graph(){
    sglib::OutputLog(sglib::INFO) << "Indexing graph, Counting..."<<std::endl;
    const int k = 31;
    uint64_t total_k=0;
    for (auto &n:sg.nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
    graph_kmers.clear();
    graph_kmers.reserve(total_k);
    FastaRecord r;
    KmerCountFactory<FastaRecord>kcf({k});
    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(graph_kmers);
        }
    }

    sglib::OutputLog(sglib::INFO)<<graph_kmers.size()<<" kmers in total"<<std::endl;
    sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
    std::sort(graph_kmers.begin(),graph_kmers.end());
    sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
    auto wi=graph_kmers.begin();
    auto ri=graph_kmers.begin();
    while (ri<graph_kmers.end()){
        if (wi.base()==ri.base()) ++ri;
        else if (*wi<*ri) {++wi; *wi=*ri;++ri;}
        else if (*wi==*ri){wi->merge(*ri);++ri;}
    }

    graph_kmers.resize(wi+1-graph_kmers.begin());
    sglib::OutputLog(sglib::INFO)<<graph_kmers.size()<<" kmers in index"<<std::endl;
    //TODO: remove kmers with more than X in count

//    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
//    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
//    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::reindex_graph(){
    const int k = 31;
    const int max_coverage = 1;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, {max_mem, 0, max_coverage,
                                                                              output_prefix});



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
    int counts=0;
    for (uint64_t i=0;i<graph_kmers.size();++i) {
        if (graph_kmers[i].count > 0){
            counts+=1;
        }
        kmer_map[graph_kmers[i].kmer] = i;
    }
    std::cout << "loaded " << graph_kmers.size() << " graph kmers " << std::to_string(counts) << " nozero counts" << std::endl;
    //read-to-node
    uint64_t ccount;
    inf.read(( char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        read_counts.emplace_back();
        read_counts.back().resize(kcount);
        inf.read(( char *) read_counts.back().data(), sizeof(uint16_t) * kcount);
    }

}

void KmerCompressionIndex::save_to_disk(std::string filename, int lib) {
    std::ofstream of(filename);
    //read-to-tag
    uint64_t kcount=graph_kmers.size();
    of.write((const char *) &kcount,sizeof(kcount));
    of.write((const char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
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
}

void KmerCompressionIndex::start_new_count(){
    read_counts.emplace_back();
    read_counts.back().resize(graph_kmers.size(),0);
}

void KmerCompressionIndex::add_counts_from_file(std::string filename) {

    std::cout<<"counting from file, graph_kmers.size(): "<< graph_kmers.size() << " sg.nodes.size() " << sg.nodes.size() <<std::endl;

    FastqReader<FastqRecord> fastqReader({0},filename);
    std::atomic<uint64_t> present(0), absent(0), rp(0);

    for (uint64_t i=0;i<graph_kmers.size();++i) kmer_map[graph_kmers[i].kmer]=i;
    size_t  sum = 0;
    size_t  kmers_in_reads = 0;
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

            for (uint64_t i=0;i<read_counts.back().size();++i) {
                read_counts.back()[i]+=thread_counts[i];
                sum += thread_counts[i];
                if (thread_counts[i] != 0){
                    kmers_in_reads += 1;
                }
            }
        }
    }

    std::cout << rp << " reads processed "<< present <<" / " << present+absent << " kmers found, mean kmers in reads: " << sum/kmers_in_reads<< std::endl;
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
    mean = ((double)tuniq)/cuniq;
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



double KmerCompressionIndex::compute_compression_for_unique_node(sgNodeID_t _node, uint16_t max_graph_freq) {
    const int k=31;
    auto n=_node>0 ? _node:-_node;
    auto & node=sg.nodes[n];

    //eliminate "overlapping" kmers
    int32_t max_bw_ovlp=0;
    int32_t max_fw_ovlp=0;
    for (auto bl:sg.get_bw_links(n)) {
        if (bl.dist<0){
            auto ovl=-bl.dist+1-k;
            if (ovl>max_bw_ovlp) max_bw_ovlp=ovl;
        }
    }
    for (auto fl:sg.get_fw_links(n)) {
        if (fl.dist<0){
            auto ovl=-fl.dist+1-k;
            if (ovl>max_fw_ovlp) max_fw_ovlp=ovl;
        }
    }
    int64_t newsize=node.sequence.size();
    newsize=newsize-max_bw_ovlp-max_fw_ovlp;
    //if (n/10==50400){
    //    std::cout<<"node "<<n<<" size="<<node.sequence.size()<<" max_bw_olv="<<max_bw_ovlp<<" max_fw_ovl="<<max_fw_ovlp<<" newlength="<<newsize<<std::endl;
    //}
    if (newsize<k) return ((double)0/0);
    auto s=node.sequence.substr(max_bw_ovlp,newsize);
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(s,k);


    skf.create_kmers(nkmers);

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (nk!=graph_kmers.end() and nk->kmer == kmer and nk->count<=max_graph_freq) {
            kcount+=nk->count;
            kcov+=read_counts[0][nk-graph_kmers.begin()];
        }
    }

    return (((double) kcov)/kcount )/uniq_mode;
}


double KmerCompressionIndex::compute_compression_for_node_old(sgNodeID_t _node, uint16_t max_graph_freq, int dataset) {

    auto & node=sg.nodes[_node>0 ? _node:-_node];
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence,31);
    skf.create_kmers(nkmers);
    //std::cout << node.sequence << std::endl;
    int counter = 0;

    uint64_t kcount=0,kcov=0;
    //std::cout << "kmers in nodes: " << nkmers.size() << std::endl;
    for (auto &kmer : nkmers){
        // find kmer in graph kmer with count > 0?
        auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (nk!=graph_kmers.end() and nk->kmer == kmer and nk-> count > 0) {
            counter +=1;
            ++kcount;// inrement number of (unique??- now removed count = 1 ) kmers on node
            kcov+=read_counts[dataset][nk-graph_kmers.begin()]; // inrement coverage by count for this kmer in read set
        }

    }
    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return (((double) kcount) )/nkmers.size();
}




double KmerCompressionIndex::compute_kcov_for_node(sgNodeID_t _node) {

    auto & node=sg.nodes[_node>0 ? _node:-_node];
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
        if (graph_kmers[kmer_map[kmer]].count > 0) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            ++kcount;// inrement number of (unique??- now removed count = 1 ) kmers on node
            kcov+=read_counts[current_lib][kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }
    //std::cout << "kcount: " << kcount << " kcov " << kcov << " kcountcount: " << kcountcount << " kcountu: " << kcountu << " kcovu " << kcovu << " kcountcountu: " << kcountcountu <<std::endl;

    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    return kcov;
}


std::vector<double> KmerCompressionIndex::compute_compression_for_node(sgNodeID_t _node, uint16_t max_graph_freq, int dataset) {

    auto & node=sg.nodes[_node>0 ? _node:-_node];
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence,31);
    skf.create_kmers(nkmers);
    //std::cout << node.sequence << std::endl;
    int counter = 0;

    uint64_t kcount=0,kcov=0,kcountcount=0,kcountu=0,kcovu=0,kcountcountu=0, counteru=0;
    //std::cout << "kmers in node: "<< _node << ", " << nkmers.size() << std::endl;
    for (auto &kmer : nkmers){
        // find kmer in graph kmer with count > 0?
        // need index of kmer in hraph_kmera - must be a better eay

        // n o idea what i was doing there... it copied from abive...
        //auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (graph_kmers[kmer_map[kmer]].count > 0) {// should scale non uniwue kmers by number occurnces in graph
            counter +=1;
            kcountcount += graph_kmers[kmer_map[kmer]].count;
            ++kcount;// inrement number of (unique??- now removed count = 1 ) kmers on node
            kcov+=read_counts[dataset][kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }
        if (graph_kmers[kmer_map[kmer]].count == 1) {
            counteru +=1;
            kcountcountu += graph_kmers[kmer_map[kmer]].count;
            ++kcountu;// inrement number of (unique??- now removed count = 1 ) kmers on node
            kcovu+=read_counts[dataset][kmer_map[kmer]]; // inrement coverage by count for this kmer in read set
        }

    }
    //std::cout << "kcount: " << kcount << " kcov " << kcov << " kcountcount: " << kcountcount << " kcountu: " << kcountu << " kcovu " << kcovu << " kcountcountu: " << kcountcountu <<std::endl;

    // number of times kmers in this node appear in reads, scaled by mod coverage of unique kmers
    std::vector<double> res = {kcount, kcov, kcountcount,kcountu, kcovu, kcountcountu, nkmers.size(), kcov/nkmers.size()};
    return res;
}