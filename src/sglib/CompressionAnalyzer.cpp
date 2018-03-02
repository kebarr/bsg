//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_name(outfile_prefix + ".txt"), outfile_csv_name(outfile_prefix+".csv"),kci(InitializeKCI()){
    InitializeKCI();
};

void CompressionAnalyzer::InitializeLib(std::string lib_name_r1, std::string lib_name_r2) {
    std::cout << "Initializing  compression analysis of " << lib_name_r1 << std::endl;

    NodeCompressions nc = {lib_name_r1, lib_name_r2};
    nc.compressions.resize(sg.nodes.size(), -1);
    compressions.push_back(nc);
    kci.start_new_count();
    kci.add_counts_from_file(lib_name_r1);
    kci.add_counts_from_file(lib_name_r2);

};

void CompressionAnalyzer::CalculateCompressions() {
    std::cout << "Calvulating compressions for each of " << compressions.size() << " read sets" << std::endl;

    for (auto nc:compressions) {
        Calculate(nc);
    }
}


std::vector<double > CompressionAnalyzer::CompressionStats(std::vector<double> res) {
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

KmerCompressionIndex CompressionAnalyzer::InitializeKCI () {
    std::cout << "Initializing kmer copression index, indexig sequene graph" << std::endl;
    KmerCompressionIndex kci(sg, max_mem_gb * 1024L * 1024L * 1024L);

    kci.index_graph();

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

// can wtite and test many versions of these- e.g taking inti acoount average compression for all contigs,
// allow for >2 repeats, vary heuristics and heristic parametes
std::vector<double> CompressionAnalyzer::AnalyseRepeat(std::vector<double> repeat_compressions, double tolerance=0.95, double diff_threshold=10) {
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
    if (exit)
        return  {-1};
    // see if repeat phased by reads- if compressions on each side pair to be reasonably close
    // sum of compressions on each side should be a multiple of middle compression

    auto in_sum = repeat_compressions[1] + repeat_compressions[2];
    auto out_sum = repeat_compressions[3] + repeat_compressions[4];
    double in_out_sane = 0;
    double repeat_count_sane = 0;

    std::vector<double> res;
    if (tolerance*in_sum < out_sum < (1-tolerance)*in_sum){
        in_out_sane = 1;
        // not sure this actually works for way i\m calculating 'compression'
        // contig repeated 5 timea should have 5*kmers in reads than average, and 5*reads going in, split in a sane way- i/e. shouldn't be 1 kmer on one in, 100 on other in, then 50/50 out
        if (std::abs(in_sum - repeat_compressions[0]) < tolerance && std::abs(out_sum - repeat_compressions[0]) < tolerance){
            repeat_count_sane = 1;
        }
        // in this ase check if resolves repeat, find out closest to in and see if close enough to call
        auto pairs = repeat_compressions[3]-repeat_compressions[1] < repeat_compressions[4]-repeat_compressions[1] ? std::make_pair(3, 4) : std::make_pair(4, 3);
        if (std::abs(repeat_compressions[1] - repeat_compressions[std::get<0>(pairs)]) < std::abs((repeat_compressions[1] - repeat_compressions[std::get<1>(pairs)]*diff_threshold)) && std::abs(repeat_compressions[2] - repeat_compressions[std::get<1>(pairs)]) < std::abs((repeat_compressions[2] - repeat_compressions[std::get<1>(pairs)]*diff_threshold))){
            // then accprding to this arbitrary heiristic, we resolve to get 0 with pair 0
            res.push_back(repeat_compressions[std::get<0>(pairs)]);
            res.push_back(repeat_compressions[1] + repeat_compressions[std::get<0>(pairs)]);

            res.push_back(repeat_compressions[std::get<1>(pairs)]);
            res.push_back(repeat_compressions[2] + repeat_compressions[std::get<1>(pairs)]);

        } else {
            // if its not resolved its more useful to know how different they were
            res.push_back(0);
            res.push_back(repeat_compressions[1] - repeat_compressions[std::get<0>(pairs)]);
            res.push_back(0);
            res.push_back(repeat_compressions[2] - repeat_compressions[std::get<1>(pairs)]);

        }
        res.push_back(in_out_sane);
        res.push_back(repeat_count_sane);
    }
    return  res;

}


void CompressionAnalyzer::Calculate(NodeCompressions & nc){
    std::ofstream outfile;
   outfile.open(outfile_name, std::ofstream::out |std::ofstream::app);
    std::ofstream outfile_csv;
    outfile_csv.open(outfile_csv_name, std::ofstream::out |std::ofstream::app);
    std::vector<std::vector<sgNodeID_t > > repeats;
    std::vector<double > repeat_contig_values;
    std::vector<double > all_repeat_contig_values;
    int count =0;
    double resolved_repeat_count =0;
    double  in_out_sane =0;
    int repeated_contig_sane =0;

    for (sgNodeID_t counter = 0; counter < sg.nodes.size(); counter++) {
                if (nc.compressions[counter] == -1) {
                    if (sg.is_canonical_repeat(counter)) {
                        std::vector<sgNodeID_t > repeat_contigs = {counter};
                        int nonzeros = 0;
                        auto bw = sg.get_bw_links(counter);
                        auto fw = sg.get_fw_links(counter);
                        // percent present/absent doesn't do it - or not obviously, now dropped unique kmer requirement
                        for (auto b: bw) {
                            auto kci_node = kci.compute_compression_for_node(b.dest, 10, nc.index);
                            outfile_csv << kci_node << ", ";
                            auto ind = b.dest > 0 ? b.dest : -b.dest;
                            outfile << sg.oldnames[ind] << ": " << kci_node << ", "
                                          << sg.nodes[ind].sequence.size() << ", ";

                            nc.compressions[b.dest] = kci_node;
                            if (kci_node > 0) nonzeros += 1;
                            repeat_contigs.push_back(b.dest);
                            all_repeat_contig_values.push_back(kci_node);
                            kci_node = 0;

                        }
                        auto kci_node = kci.compute_compression_for_node(counter, 10, nc.index);
                        outfile_csv << kci_node << ", ";
                        repeat_contig_values.push_back(kci_node);
                        if (kci_node > 0) nonzeros += 1;
                        nc.compressions[counter] = kci_node;
                        all_repeat_contig_values.push_back(kci_node);

                        kci_node = 0;

                        outfile << sg.oldnames[counter] << ": " << kci_node << ", "
                                      << sg.nodes[counter].sequence.size() << ", ";
                        for (auto f: fw) {
                            auto kci_node = kci.compute_compression_for_node(f.dest, 10, nc.index);
                            outfile_csv << kci_node << ", ";
                            auto ind = f.dest > 0 ? f.dest : -f.dest;
                            outfile << sg.oldnames[ind] << ": " << kci_node << ", "
                                          << sg.nodes[ind].sequence.size() << ", ";
                            if (kci_node > 0) nonzeros += 1;


                            nc.compressions[f.dest] = kci_node;
                            repeat_contigs.push_back(f.dest);
                            all_repeat_contig_values.push_back(kci_node);
                            repeat_contig_values.push_back(kci_node);

                            kci_node = 0;

                        }
                        outfile << std::endl;
                        count += 1;
                        if (nonzeros >= 3) {
                            auto res = AnalyseRepeat(all_repeat_contig_values);
                            in_out_sane += res[4];
                            repeated_contig_sane += res[5];

                            if (res[0] != 0) {
                                outfile << "Resolved: ";
                                resolved_repeat_count += 1;
                            }

                            for (int r= 0;  r <res.size() ; r++) {
                                outfile << r << ", ";
                            }
                            //
                            nc.canonical_repeats.push_back(repeat_contigs);
                        } else {
                            outfile << "repeat not present " << std::endl;
                        }
                        outfile << std::endl;

                        all_repeat_contig_values.clear();
                    } else {
                        auto kci_node = kci.compute_compression_for_node(counter, 10, nc.index);
                        outfile_csv << kci_node << ", ";
                        nc.compressions[counter] = kci_node;
                    }
                }
    }

outfile_csv << std::endl;
    auto all_stats = CompressionStats(nc.compressions);//  {sum, sd, mean, *max, *min};

std::cout << "stats for all contigs, min  " << all_stats[4]<< " max: " << all_stats[3] << " mean: "<< all_stats[2]
<<   " stdev: " << all_stats[1]<< std::endl;
    std::vector<double > repeat_compressions;
    auto repeat_stats = CompressionStats(repeat_contig_values);

    std::cout << "stats for all contigs in cononical repeats, min  " << repeat_stats[4]<< " max: " << repeat_stats[3] << " mean: "<< repeat_stats[2]
              <<   " stdev: " << repeat_stats[1]<< std::endl;


    auto all_repeat_stats = CompressionStats(all_repeat_contig_values);


    std::cout << "stats for all contigs in cononical repeats, min  " << all_repeat_stats[4]<< " max: " << all_repeat_stats[3] << " mean: "<< all_repeat_stats[2]
              <<   " stdev: " << all_repeat_stats[1]<< std::endl;

    std::vector<double > all_repeat_compressions;
    for (auto r:nc.canonical_repeats) {
        // acytually repeated node always first
        repeat_compressions.push_back(nc.compressions[r[0]]);
        for(auto c: r) {
            all_repeat_compressions.push_back(nc.compressions[c]);
        }
    }

    std::cout << "calculated compression for " << nc.lib_name_r1 << " for " <<
              compressions.size() << " nodes, maps to " << repeat_contig_values.size() << " repeats, of which " << resolved_repeat_count << "resolved. \n"
                      " calxulated scores for  " << repeat_contig_values.size() << " repeated contigs, including the in/out contigs,  " << all_repeat_contig_values.size() << " scored in total, of which  " << in_out_sane << " in/out compression sums are consistent and " << repeated_contig_sane << " repeated contig compressions are consistent with in/out suns"  << std::endl;


        }
