//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_name(outfile_prefix + ".txt"), outfile_csv_name(outfile_prefix+".csv"),kci(this->sg, this->max_mem_gb*1024L*1024L*1024L){
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

void CompressionAnalyzer::CalculateCompressions() {
    std::cout << "Calvulating compressions for each of " << compressions.size() << " read sets" << std::endl;

    for (auto &nc:compressions) {
        Calculate(nc);
    }

    for (auto &nc:compressions) {
        std::cout << nc.index << ": ";
        for (auto l:nc.compressions) {
            std::cout << l << ", ";
        }
        std::cout << std::endl;
        for (auto s:nc.canonical_repeats) {
            for (auto l:s){
            std::cout << sg.oldnames[l > 0? l: -l] << ", ";
            }
        }
        std::cout << std::endl;

    }


}



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

void CompressionAnalyzer::InitializeKCI () {
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


void CompressionAnalyzer::Calculate(NodeCompressions & nc){
    std::cout << "calculating compressions of " << nc.lib_name_r1 << " and " << nc.lib_name_r2 << " writing to " << outfile_name<< std::endl;
    std::ofstream outfile;
   outfile.open(outfile_name, std::ofstream::out |std::ofstream::app);
    std::ofstream outfile_csv,outfile_csv1, outfile_csv2,outfile_csv3,outfile_csv4,outfile_csv5,outfile_csv6,outfile_csv7;
    outfile_csv.open(outfile_csv_name, std::ofstream::out |std::ofstream::app);
    outfile_csv1.open(outfile_csv_name1, std::ofstream::out |std::ofstream::app);
    outfile_csv2.open(outfile_csv_name2, std::ofstream::out |std::ofstream::app);
    outfile_csv3.open(outfile_csv_name3, std::ofstream::out |std::ofstream::app);
    outfile_csv4.open(outfile_csv_name4, std::ofstream::out |std::ofstream::app);
    outfile_csv5.open(outfile_csv_name5, std::ofstream::out |std::ofstream::app);
    outfile_csv6.open(outfile_csv_name6, std::ofstream::out |std::ofstream::app);
    outfile_csv7.open(outfile_csv_name7, std::ofstream::out |std::ofstream::app);

    std::vector<std::ofstream> csvs;
        csvs.push_back(std::move(outfile_csv));
    csvs.push_back(std::move(outfile_csv1));

    csvs.push_back(std::move(outfile_csv2));
    csvs.push_back(std::move(outfile_csv3));
    csvs.push_back(std::move(outfile_csv4));
    csvs.push_back(std::move(outfile_csv5));
    csvs.push_back(std::move(outfile_csv6));
    csvs.push_back(std::move(outfile_csv7));

   // = {outfile_csv, outfile_csv1,outfile_csv2, outfile_csv3, outfile_csv4, outfile_csv5, outfile_csv6};

    std::vector<std::vector<double > >repeat_contig_values;
    std::vector<std::vector<double > > all_repeat_contig_values;
    std::vector<sgNodeID_t > repeat_contigs = {};
    std::vector< std::vector<sgNodeID_t > > canonical;
    int count =0;
    int repeats = 0;
    std::map<int, std::string> metrics = {{0, "kcount"}, {1,"kcov"} ,{2, "kcountcount"}, {3,"kcountu"}, {4, "kcovu"}, {5,"kcountcount"}, { 6, "nkmers.size()"}, {7, "kcov/nkmers.size()"}};
std::vector<sgNodeID_t > resolved_repeat_indices;
    std::map<int, int> resolved_repeat_count ={{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0}};
    std::map<int, int>  in_out_sane ={{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0}};
    std::map<int, int> repeated_contig_sane ={{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0}};

    for (sgNodeID_t counter = 1; counter < sg.nodes.size(); counter++) {
        /*std::cout << "Counter: " << counter  << std::endl;
        std::cout << "nc.compressionssize " << nc.compressions.size() << std::endl;
        if (counter < sg.oldnames.size()) {

            std::cout << " sg.oldnames: " << sg.oldnames[counter] << std::endl;
        }
        std::cout << " nc " << nc.lib_name_r2 << std::endl;
        std::cout << " r1: " << nc.lib_name_r1 << std::endl;
        std::cout  << " kci.read_counts.size() "<< kci.read_counts.size()<< std::endl;
        std::cout << " ind: "<< nc.index << std::endl; std::cout << " nc.canonical_repeats.si" << nc.canonical_repeats.size() << std::endl;
*/
        // lines4 pnted "kmers in node " 33 times so get info about roughly where that is
                    //std::cout << "Counter: " << counter << " sg.oldnames: " << sg.oldnames[counter] << " nc " << nc.lib_name_r2 << " r1: " << nc.lib_name_r1 << " kci.read_counts.size() "<< kci.read_counts.size() << " ind: "<< nc.index << " nc.canonical_repeats.si" << nc.canonical_repeats.size() << std::endl;


        std::vector<std::vector<double > >local_repeat_contig_values;

        auto kci_node = kci.compute_compression_for_node(counter, 10, nc.index);
        int nonzeros = 0;
        if (kci_node[0] > 0) nonzeros += 1;
        all_repeat_contig_values.push_back(kci_node);
        int c = 0;// always record kci in csv
        for (auto k: kci_node) {
            csvs[c] << k << ", ";
            c+=1;
        }
        local_repeat_contig_values.push_back(kci_node);


        if (nc.compressions[counter] == -1) {

                    if (sg.is_canonical_repeat(counter)) {
                        count += 1;

                        repeats += 1;
                        std::vector<sgNodeID_t > repeat_contigs = {counter};

                        auto bw = sg.get_bw_links(counter);
                        auto fw = sg.get_fw_links(counter);

                        for (auto b: bw) {
                            auto kci_node = kci.compute_compression_for_node(b.dest, 10, nc.index);
                            count += 1;

                            auto ind = b.dest > 0 ? b.dest : -b.dest;
                            std::cout << "node: " << sg.oldnames[ind] << " ";
int c=0;
                            for (auto k: kci_node) {
                                csvs[c] << k << ", ";
                                c+=1;
                            }

                            nc.compressions[b.dest] = kci_node[1]/kci_node[6];
                            if (kci_node[0] > 0) nonzeros += 1;
                            repeat_contigs.push_back(ind);
                            local_repeat_contig_values.push_back(kci_node);
                            all_repeat_contig_values.push_back(kci_node);

                        }





                        for (auto f: fw) {
                            auto kci_node = kci.compute_compression_for_node(f.dest, 10, nc.index);
                            count += 1;
                            auto ind = f.dest > 0 ? f.dest : -f.dest;
int c=0;
                            for (auto k: kci_node) {
                                csvs[c] << k << ", ";

                                c+=1;
                            }
                            if (kci_node[0] > 0) nonzeros += 1;


                            nc.compressions[ind] =  kci_node[1]/kci_node[6];
                            repeat_contigs.push_back(ind);
                            all_repeat_contig_values.push_back(kci_node);
                            local_repeat_contig_values.push_back(kci_node);



                        }
                        outfile << std::endl;
                        if (nonzeros >= 3) {
                            auto res = AnalyseRepeat(local_repeat_contig_values);
for (int i=0; i < res.size() ; i++) {// i is index of node

    in_out_sane[i] += res[i][4];
    repeated_contig_sane[i] += res[i][5];
    if (res[i][0] != 0) {
        outfile << i << ": ";
        std::cout << i << ": ";
        resolved_repeat_indices.push_back(counter);

        outfile << "\nResolved: " << i << "  " << metrics[i] << " ";
        std::cout << "\nResolved: " << i << " " << metrics[i] << " ";


        outfile << std::endl;
        std::cout << std::endl;
        resolved_repeat_count[i] += 1;

        for (int r = 0; r < res[i].size(); r++) {
            outfile << metrics[r] << ": " << res[i][r] << ", ";
            std::cout << metrics[r] << ": " << res[i][r] << ", ";

        }


        outfile << std::endl;
        std::cout << std::endl;
    }
}
                            std::cout << "repeats";
                            for (auto r: repeat_contigs){
                                if (r < sg.oldnames.size()) {
                                    std::cout << r << ": "  << sg.oldnames[r] << ", ";

                                } else {
                                    std::cout << r << ": no old name"  << ", ";

                                }
                            }

                            nc.canonical_repeats.push_back({});

                            for (auto r: repeat_contigs){
                                nc.canonical_repeats.back().push_back(r);
                            }
                                      //
                            canonical.push_back(repeat_contigs);
                        } else {
                            outfile << "repeat not present " << std::endl;
                        }
                        outfile << std::endl;


                    }
                }
        nc.compressions[counter] = kci_node[1]/kci_node[6];

        local_repeat_contig_values.clear();
    }
    outfile_csv << std::endl;
    outfile_csv1 << std::endl;
    outfile_csv2 << std::endl;
    outfile_csv3 << std::endl;
    outfile_csv4 << std::endl;
    outfile_csv5 << std::endl;
    outfile_csv6 << std::endl;

    auto all_stats = CompressionStats(nc.compressions);//  {sum, sd, mean, *max, *min};

std::cout << "stats for all contigs, min  " << all_stats[4]<< " max: " << all_stats[3] << " mean: "<< all_stats[2]
<<   " stdev: " << all_stats[1]<< std::endl;
    std::vector<double > repeat_compressions;
    /*
    for (int i=0; i < repeat_contig_values.size() ; i++) {
        auto repeat_stats = CompressionStats(repeat_contig_values[i]);

        std::cout << "stats for all contigs in cononical repeats, min  " << repeat_stats[4] << " max: "
                  << repeat_stats[3] << " mean: " << repeat_stats[2]
                  << " stdev: " << repeat_stats[1] << std::endl;


        auto all_repeat_stats = CompressionStats(all_repeat_contig_values[i]);


        std::cout << "stats for all contigs in cononical repeats, min  " << all_repeat_stats[4] << " max: "
                  << all_repeat_stats[3] << " mean: " << all_repeat_stats[2]
                  << " stdev: " << all_repeat_stats[1] << std::endl;
    }*/
    std::vector<double > all_repeat_compressions;
    for (auto r:nc.canonical_repeats) {
        // acytually repeated node always first
        repeat_compressions.push_back(nc.compressions[r[0]]);
        std::cout << "repeated: " << sg.oldnames[r[0]];
        for(auto c: r) {
            all_repeat_compressions.push_back(nc.compressions[c]);
            std::cout << sg.oldnames[c] << " " << nc.compressions[c] ;


        }
        std::cout <<std::endl;
    }


    int total_resolved = 0;
    int best_metric = -1;
    int most_resolved = -1;
    for (auto r:resolved_repeat_count){
        if (r.second > most_resolved){
            most_resolved = r.second;
            best_metric = r.first;
        }
        if (r.second != 0){}

    }
    std::cout << "best_metric: " << metrics[best_metric] << " resolves " << most_resolved << std::endl;
    std::cout << "calculated compression for " << nc.lib_name_r1 << " for " <<
              compressions.size() << " nodes, maps to more than 3 nodes of " << nc.canonical_repeats.size() << " repeats, of which " << resolved_repeat_indices.size() << "resolved. \n"<<
                     " including the in/out contigs,  " << all_repeat_contig_values.size() << " scored in total, of which  " << in_out_sane[best_metric] << " in/out compression sums are consistent and " << repeated_contig_sane[best_metric] << " repeated contig compressions are consistent with in/out suns"  << std::endl;


    return;

    //TODO:: do this properly!! make tedt datast with canibical repeats + read sets resolving!!!
        }



