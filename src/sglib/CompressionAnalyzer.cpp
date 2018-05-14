//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"
#include "RepeatAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_name(outfile_prefix + ".txt"), outfile_csv_name(outfile_prefix+".csv"),kci(this->sg, this->max_mem_gb*1024L*1024L*1024L), ra(this->sg){
    std::vector<std::string> csvs = {outfile_csv_name1,outfile_csv_name2,outfile_csv_name3,outfile_csv_name4,outfile_csv_name5,outfile_csv_name6,outfile_csv_name7};
    for(int i=1; i <7; i++){
        csvs[i] = outfile_prefix+std::to_string(i) + ".csv";
    }
    this->outfile_csv_name1 = csvs[0];
    this->outfile_csv_name2 = csvs[1];
    this->outfile_csv_name3 = csvs[2];
    this->outfile_csv_name4 = csvs[3];
    this->outfile_csv_name5 = csvs[4];
    this->outfile_csv_name6 = csvs[5];
    this->outfile_csv_name7 = csvs[6];
    InitializeKCI();
    ra.FindRepeats();

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



std::vector<double > CompressionAnalyzer::CompressionStats(std::vector<double> res) {
    if (res.size() == 0 ){
        return {0,0,0,0,0};
    }
    auto sum = std::accumulate(res.begin(), res.end(), 0);
    double sd;
    double mean = sum / res.size();
    for (auto c:res) {
        sd += (c - mean) * (c - mean);
    }
    sd = (std::pow(sd, 0.5)) / ( res.size() - 1);
    auto max = std::max_element(res.begin(), res.end());
    auto min = std::min_element(res.begin(), res.end());
    std::vector<double > final = {sum, sd, mean, *max, *min};
    return final;

};


std::vector<Repeat> CompressionAnalyzer::FindGraphRepeats(){
    /*
}
    std::vector<Repeat> repeats;
    size_t max_deg = 0;
    size_t  min_deg = sg.nodes.size();
    for (sgNodeID_t counter = 1; counter < sg.nodes.size(); counter++) {

 if (sg.is_canonical_repeat(counter)) {
    Repeat rep;
     rep.repeated_contig = counter;
     auto bw = sg.get_bw_links(counter);
     auto fw = sg.get_fw_links(counter);
     for (auto b:bw) rep.in_contigs.push_back(b.dest);
     for (auto f:fw) rep.in_contigs.push_back(f.dest);
     rep.degree = bw.size();
     if (bw.size() < min_deg) min_deg = bw.size();

     if (bw.size() > max_deg) max_deg = bw.size();
    repeats.push_back(rep);
 }
        std::cout << "counted " << repeats.size() << " repeats with min degree  " << min_deg << " and max degree " << max_deg <<std::endl;
        return  repeats;*/
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
    std::ofstream outfile_csv1;

    std::ofstream outfile_csv2;

    std::ofstream outfile_csv3;
    std::ofstream outfile_csv4;
    std::ofstream outfile_csv5;
    std::ofstream outfile_csv6;
    std::ofstream outfile_csv7;

    outfile_csv.open(outfile_csv_name, std::ofstream::out );
    outfile_csv1.open(outfile_csv_name1, std::ofstream::out );
    outfile_csv2.open(outfile_csv_name2, std::ofstream::out );
    outfile_csv3.open(outfile_csv_name3, std::ofstream::out );
    outfile_csv4.open(outfile_csv_name4, std::ofstream::out );
    outfile_csv5.open(outfile_csv_name5, std::ofstream::out );
    outfile_csv6.open(outfile_csv_name6, std::ofstream::out );
    outfile_csv7.open(outfile_csv_name7, std::ofstream::out );

    for (size_t counter = 0; counter < sg.nodes.size(); counter++) {
        outfile_csv << sg.oldnames[counter] << ", ";
        outfile_csv1 << sg.oldnames[counter] << ", ";
        outfile_csv2 << sg.oldnames[counter] << ", ";
        outfile_csv3 << sg.oldnames[counter] << ", ";
        outfile_csv4 << sg.oldnames[counter] << ", ";
        outfile_csv5 << sg.oldnames[counter] << ", ";
        outfile_csv6<< sg.oldnames[counter] << ", ";

        outfile_csv7<< sg.oldnames[counter] << ", ";

    }
    outfile_csv << std::endl;
    outfile_csv1 << std::endl;
    outfile_csv2 << std::endl;
    outfile_csv3 << std::endl;
    outfile_csv4 << std::endl;
    outfile_csv5 << std::endl;
    outfile_csv6 << std::endl;
    outfile_csv7 << std::endl;

    std::cout << "outfile_csv_name6 " << outfile_csv_name6 << std::endl;
    if (kci.read_counts.size()>0) {
        kci.compute_compression_stats();
        kci.dump_histogram(outfile_name + "kci_histogram.csv");
    }

}

// can wtite and test many versions of these- e.g taking inti acoount average compression for all contigs,
// allow for >2 repeats, vary heuristics and heristic parametes
std::vector<std::vector<double>> CompressionAnalyzer::AnalyseRepeat(std::vector<std::vector<double>> repeat_compressions, double tolerance=0.8, double diff_threshold=0.8) {

  bool exit = false;

    if (repeat_compressions.size() > 5) {
        std::cout << "repeats of more thsan 2 not yet supported " << std::endl;

        if (repeat_compressions.size() % 2 != 1) {

            std::cout << " even sized vector of  " << repeat_compressions.size()
                      << " contig ids is ont a valid representation of a canonical repeat" << std::endl;
        }
        exit = true;
    } else if (repeat_compressions.size() < 5) {
        std::cout << "vector of size " << repeat_compressions.size() << " thsan 2 not yet supported " << std::endl;
        if (repeat_compressions.size() % 2 != 1) {

            std::cout << " even sized vector of  " << repeat_compressions.size()
                      << " contig is is ont a valid representation of a canonical repeat" << std::endl;
        }
        exit = true;
    }
    if (exit){
        return  {{-1}};
    }
    std::vector<std::vector<double>> res_all;
    // test each metric
    for (int i=0; i <repeat_compressions[0].size(); i++) {

        bool exit = false;
        // see if repeat phased by reads- if compressions on each side pair to be reasonably close
        // sum of compressions on each side should be a multiple of middle compression

        auto in_sum = repeat_compressions[1][i] + repeat_compressions[2][i];
        auto out_sum = repeat_compressions[3][i] + repeat_compressions[4][i];
        double in_out_sane = 0;
        double repeat_count_sane = 0;

        std::vector<double> res = {0, 0, 0, 0, 0, 0};
        double ratioinout = in_sum< out_sum ? in_sum/out_sum : out_sum/in_sum;
        std::cout << "in: " << repeat_compressions[1][i] << " " << repeat_compressions[2][i] << " in sum: " << in_sum << " out: " << repeat_compressions[3][i] << repeat_compressions[4][i] << " out sum: " << out_sum <<  " ratioinout " << ratioinout <<std::endl;


        if (ratioinout > tolerance) {
            in_out_sane = 1;
            auto a = std::abs(out_sum - repeat_compressions[0][i]) < tolerance;
            auto b = std::abs(in_sum - repeat_compressions[0][i]) > tolerance;

            // not sure this actually works for way i\m calculating 'compression'
            double ratio1 = in_sum < repeat_compressions[0][i] ? in_sum / repeat_compressions[0][i] :
                            repeat_compressions[0][i] / in_sum;
            double ratio2 = out_sum < repeat_compressions[0][i] ? out_sum / repeat_compressions[0][i] :
                            repeat_compressions[0][i] / out_sum;

            // contig repeated 5 timea should have 5*kmers in reads than average, and 5*reads going in, split in a sane way- i/e. shouldn't be 1 kmer on one in, 100 on other in, then 50/50 out
            if (ratio1 > tolerance && ratio2 > tolerance) {
                repeat_count_sane = 1;
            }
            // in this ase check if resolves repeat, find out closest to in and see if close enough to call
            auto pairs =
                    repeat_compressions[3][i] - repeat_compressions[1][i] <
                    repeat_compressions[4][i] - repeat_compressions[1][i]
                    ? std::make_pair(3, 4) : std::make_pair(4, 3);
            auto ratio_in_out_match1 = repeat_compressions[1][i] > repeat_compressions[std::get<0>(pairs)][i] ?
                                       repeat_compressions[std::get<0>(pairs)][i] / repeat_compressions[1][i] :
                                       repeat_compressions[1][i] / repeat_compressions[std::get<0>(pairs)][i];

            auto ratio_in_out_match2 = repeat_compressions[2][i] > repeat_compressions[std::get<1>(pairs)][i] ?
                                       repeat_compressions[std::get<1>(pairs)][i] / repeat_compressions[2][i] :
                                       repeat_compressions[2][i] / repeat_compressions[std::get<1>(pairs)][i];

            if (ratio_in_out_match1 > diff_threshold && ratio_in_out_match2 > diff_threshold) {
                // then accprding to this arbitrary heiristic, we resolve to get 0 with pair 0
                res[0] = repeat_compressions[std::get<0>(
                        pairs)][i];// compression of closest out contig to first in contig
                res[1] = repeat_compressions[1][i] +
                         repeat_compressions[std::get<0>(pairs)][i];// compression of both in contigs

                res[2] = repeat_compressions[std::get<1>(pairs)][i];
                res[3] = repeat_compressions[2][i] + repeat_compressions[std::get<1>(pairs)][i];

            } else {
                // if its not resolved its more useful to know how different they were
                res[1] = repeat_compressions[1][i] - repeat_compressions[std::get<0>(pairs)][i];
                res[3] = repeat_compressions[2][i] - repeat_compressions[std::get<1>(pairs)][i];

            }
            res[4] = in_out_sane;
            res[5] = repeat_count_sane;

            std::cout << "res" << i << ": ";
            for (auto r:res) {
                std::cout << r << " ";
            }

            std::cout << std::endl;
        }
        res_all.push_back(res);
    }
    return  res_all;

}


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

void CompressionAnalyzer::CalculateRepeatCompressions() {
   int resolved_repeat_indices = 0;
for (int i=0;i <kci.read_counts.size(); i++) {
    kci.current_lib = i;
   auto c= ra.compressions_for_read_set(compute_kcov_for_node);
    this->compressions[i].compressions = c;
}
// now have all compressions for each lib (though not sure its needed), and each repeat has the compressions from all libs so that we can see whether the different varieties resolve it

};
