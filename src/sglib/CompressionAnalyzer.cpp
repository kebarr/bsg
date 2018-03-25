//
// Created by Katie Barr (EI) on 02/03/2018.
//

#include "CompressionAnalyzer.h"


CompressionAnalyzer::CompressionAnalyzer(SequenceGraph &_sg, uint64_t max_mem_gb, std::string outfile_prefix) :sg(_sg), max_mem_gb(max_mem_gb), outfile_name(outfile_prefix + ".txt"), outfile_csv_name(outfile_prefix+".csv"),kci(this->sg, this->max_mem_gb*1024L*1024L*1024L){
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
std::vector<std::vector<double>> CompressionAnalyzer::AnalyseRepeat(std::vector<std::vector<double>> repeat_compressions, double tolerance=0.95, double diff_threshold=10) {
    for (auto r:repeat_compressions){
        for (auto a:r) {
            std::cout << a << " ";
        }
        std::cout << std::endl;
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
        return  {{-1}};
    }
    std::vector<std::vector<double>> res_all;
    for (int i=0; i <repeat_compressions.size(); i++) {

        bool exit = false;
        // see if repeat phased by reads- if compressions on each side pair to be reasonably close
        // sum of compressions on each side should be a multiple of middle compression

        auto in_sum = repeat_compressions[i][1] + repeat_compressions[i][2];
        auto out_sum = repeat_compressions[i][3] + repeat_compressions[i][4];
        double in_out_sane = 0;
        double repeat_count_sane = 0;

        std::vector<double> res = {0, 0, 0, 0, 0, 0};
        std::cout << "in sum: " << in_sum << " out sum: " << out_sum << std::endl;
        double ratioinout = in_sum< out_sum ? in_sum/out_sum : out_sum/in_sum;
        std::cout << "in sum: " << in_sum << " out sum: " << out_sum <<  " ratioinout" << ratioinout <<std::endl;


        if (ratioinout > tolerance) {
            in_out_sane = 1;
            auto a = std::abs(out_sum - repeat_compressions[0][i]) < tolerance;
            auto b = std::abs(in_sum - repeat_compressions[0][i]) > tolerance;
            std::cout << "i: " << i << " std::abs(in_sum - repeat_compressions[0]) " << std::abs(in_sum - repeat_compressions[0][i])
                      << "> tolerance  " << b << " std::abs(out_sum - repeat_compressions[0])  "
                      << std::abs(out_sum - repeat_compressions[0][i]) << "  < tolerance " <<
                      a << std::endl;
            // not sure this actually works for way i\m calculating 'compression'
            double ratio1 = in_sum< repeat_compressions[0][i] ? in_sum/repeat_compressions[0][i] : repeat_compressions[0][i]/in_sum;
            double ratio2 = out_sum< repeat_compressions[0][i]? out_sum/repeat_compressions[0][i] : repeat_compressions[0][i]/out_sum;
            std::cout << "ratio1: " << ratio1 << " ratio2: " << ratio2 << "\n std::abs(in_sum - repeat_compressions[0]) " << std::abs(in_sum - repeat_compressions[i][0])
                      << "> tolerance  " << b << " std::abs(out_sum - repeat_compressions[0])  "
                      << std::abs(out_sum - repeat_compressions[0][i]) << "  < tolerance " <<
                      a << std::endl;
            // contig repeated 5 timea should have 5*kmers in reads than average, and 5*reads going in, split in a sane way- i/e. shouldn't be 1 kmer on one in, 100 on other in, then 50/50 out
            if (ratio1< tolerance && ratio2 < tolerance) {
                repeat_count_sane = 1;
            }
            // in this ase check if resolves repeat, find out closest to in and see if close enough to call
            auto pairs =
                    repeat_compressions[i][3] - repeat_compressions[i][1] < repeat_compressions[i][4] - repeat_compressions[i][1]
                    ? std::make_pair(3, 4) : std::make_pair(4, 3);
            if (std::abs(repeat_compressions[i][1] - repeat_compressions[i][std::get<0>(pairs)]) <
                std::abs((repeat_compressions[i][1] - repeat_compressions[i][std::get<1>(pairs)] * diff_threshold)) &&
                std::abs(repeat_compressions[i][2] - repeat_compressions[i][std::get<1>(pairs)]) <
                std::abs((repeat_compressions[i][2] - repeat_compressions[i][std::get<1>(pairs)] * diff_threshold))) {
                // then accprding to this arbitrary heiristic, we resolve to get 0 with pair 0
                res[0] = repeat_compressions[i][std::get<0>(pairs)];// compression of closest out contig to first in contig
                res[1] = repeat_compressions[i][1] +
                         repeat_compressions[i][std::get<0>(pairs)];// compression of both in contigs

                res[2] = repeat_compressions[i][std::get<1>(pairs)];
                res[3] = repeat_compressions[i][2] + repeat_compressions[i][std::get<1>(pairs)];

            } else {
                // if its not resolved its more useful to know how different they were
                res[1] = repeat_compressions[i][1] - repeat_compressions[i][std::get<0>(pairs)];
                res[3] = repeat_compressions[i][2] - repeat_compressions[i][std::get<1>(pairs)];

            }
            res[4] = in_out_sane;
            res[5] = repeat_count_sane;
        }
        std::cout << "res: ";
        for (auto r:res) {
            std::cout << r << " ";
        }
        std::cout << std::endl;
    }
    return  res_all;

}


void CompressionAnalyzer::Calculate(NodeCompressions & nc){
    std::cout << "calculating compressions of " << nc.lib_name_r1 << " and " << nc.lib_name_r2 << std::endl;
    std::ofstream outfile;
   outfile.open(outfile_name, std::ofstream::out |std::ofstream::app);
    std::ofstream outfile_csv;
    outfile_csv.open(outfile_csv_name, std::ofstream::out |std::ofstream::app);
    std::vector<std::vector<double > >repeat_contig_values;
    std::vector<std::vector<double > > all_repeat_contig_values;
    std::vector<sgNodeID_t > repeat_contigs = {};
    std::vector< std::vector<sgNodeID_t > > canonical;
    int count =0;
    int repeats = 0;
    double resolved_repeat_count =0;
    double  in_out_sane =0;
    int repeated_contig_sane =0;
    for (auto i= 1 ; i < sg.links.size() ; i++){
        if (i < sg.oldnames.size()){
            std::cout << sg.links[i].size() << ": " << sg.oldnames[i] << " ";
        } else {
            std::cout << sg.links[i].size() << " i: " << i;
        }
        std::cout << " fw: " << sg.get_fw_links(i).size() << " bw: " << sg.get_bw_links(i).size() << std::endl;
        for (auto n:sg.links[i]) {
            if (n.source < sg.oldnames.size()){
            std::cout     << sg.oldnames[n.source > 0 ? n.source : - n.source ];
            } else {
                std::cout  << n.source    <<     "no aource ";
            }
            if (n.dest < sg.oldnames.size()){
                std::cout     << " d: " << sg.oldnames[n.dest];
            } else {
                std::cout   << n.dest  <<     "no dest ";
            }
            std::cout     << ", ";
        }
        std::cout     << std::endl;
    }
    std::cout << "nc.compressionssize " << nc.compressions.size() << std::endl;
    for (sgNodeID_t counter = 1; counter < sg.nodes.size(); counter++) {
        std::cout << "Counter: " << counter  << std::endl;
        std::cout << "nc.compressionssize " << nc.compressions.size() << std::endl;
        if (counter < sg.oldnames.size()) {

            std::cout << " sg.oldnames: " << sg.oldnames[counter] << std::endl;
        }
        std::cout << " nc " << nc.lib_name_r2 << std::endl;
        std::cout << " r1: " << nc.lib_name_r1 << std::endl;
        std::cout  << " kci.read_counts.size() "<< kci.read_counts.size()<< std::endl;
        std::cout << " ind: "<< nc.index << std::endl; std::cout << " nc.canonical_repeats.si" << nc.canonical_repeats.size() << std::endl;

        // lines4 pnted "kmers in node " 33 times so get info about roughly where that is
                    std::cout << "Counter: " << counter << " sg.oldnames: " << sg.oldnames[counter] << " nc " << nc.lib_name_r2 << " r1: " << nc.lib_name_r1 << " kci.read_counts.size() "<< kci.read_counts.size() << " ind: "<< nc.index << " nc.canonical_repeats.si" << nc.canonical_repeats.size() << std::endl;
                    for (auto e: sg.get_bw_links(counter)){
                        auto ind = e.dest > 0 ? e.dest : -e.dest;
                        std::cout << "bw " << e << " old: " << sg.oldnames[ind] <<  " comp e.dest " << nc.compressions[e.dest] << std::endl;
                        for (auto a: kci.compute_compression_for_node(e.dest, 10, nc.index) ){
                            std::cout << a << ",  ";
                        }
                        std::cout << std::endl;
                    }

                    for (auto e: sg.get_fw_links(counter)){
                        auto ind = e.dest > 0 ? e.dest : -e.dest;
                        std::cout << "fw " << e << " old: " << sg.oldnames[ind] <<   " comp e.dest " << nc.compressions[e.dest] << std::endl;
                        for (auto a: kci.compute_compression_for_node(e.dest, 10, nc.index)){
                            std::cout << a << ",  ";
                        }
                        std::cout << std::endl;
                    }


                if (nc.compressions[counter] == -1) {
                    if (sg.is_canonical_repeat(counter)) {
                        repeats += 1;
                        std::vector<sgNodeID_t > repeat_contigs = {counter};
                        std::vector<std::vector<double > >local_repeat_contig_values;
                        int nonzeros = 0;
                        auto bw = sg.get_bw_links(counter);
                        auto fw = sg.get_fw_links(counter);
                        for (auto b: bw) {
                            auto kci_node = kci.compute_compression_for_node(b.dest, 10, nc.index);
                            count += 1;

                            auto ind = b.dest > 0 ? b.dest : -b.dest;
                            for (auto k: kci_node) {
                                outfile_csv << k << ", ";
                                if (ind < sg.oldnames.size()) {
                                    outfile << sg.oldnames[ind] << ": " << k << ", "
                                            << sg.nodes[ind].sequence.size() << ", ";
                                } else {
                                    outfile << ind << " no old name: " << k << ", "
                                            << sg.nodes[ind].sequence.size() << ", ";
                                }
                            }

                            nc.compressions[b.dest] = kci_node[1]/kci_node[6];
                            if (kci_node[0] > 0) nonzeros += 1;
                            repeat_contigs.push_back(ind);
                            local_repeat_contig_values.push_back(kci_node);
                            all_repeat_contig_values.push_back(kci_node);

                        }
                        auto kci_node = kci.compute_compression_for_node(counter, 10, nc.index);
                        count += 1;
                        if (kci_node[0] > 0) nonzeros += 1;
                        nc.compressions[counter] = kci_node[1]/kci_node[6];
                        all_repeat_contig_values.push_back(kci_node);
                        local_repeat_contig_values.push_back(kci_node);

                        for (auto k: kci_node) {
                            outfile_csv << k << ", ";
                            if (counter < sg.oldnames.size()) {
                                outfile << sg.oldnames[counter] << ": " << k<< ", "
                                        << sg.nodes[counter].sequence.size() << ", ";
                            } else {
                                outfile << counter << " no old name: " << k<< ", "
                                        << sg.nodes[counter].sequence.size() << ", ";
                            }
                        }




                        for (auto f: fw) {
                            auto kci_node = kci.compute_compression_for_node(f.dest, 10, nc.index);
                            count += 1;
                            auto ind = f.dest > 0 ? f.dest : -f.dest;
                            for (auto k: kci_node) {
                                outfile_csv << k << ", ";
                                if (ind < sg.oldnames.size()) {
                                    outfile << sg.oldnames[ind] << ": " << k << ", "
                                            << sg.nodes[ind].sequence.size() << ", ";
                                } else {
                                    outfile << ind << " no old name: " << k << ", "
                                            << sg.nodes[ind].sequence.size() << ", ";
                                }
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
for (int i=0; i < res.size() ; i++) {
    std::cout << i << " ";
    in_out_sane += res[i][4];
    repeated_contig_sane += res[i][5];

    if (res[i][0] != 0) {
        outfile << "Resolved: ";
        std::cout << "Resolved: ";

        resolved_repeat_count += 1;
    }

    for (int r = 0; r < res[i].size(); r++) {
        outfile << r << ", ";
        std::cout << r << ", ";

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

                            std::cout << std::endl << "  nc.canonical_repeats " <<  nc.canonical_repeats.size() << std::endl ;
                            nc.canonical_repeats.push_back({});
                            //std::cout << "  nc.canonical_repeats 2 " <<  nc.canonical_repeats.size() << std::endl ;

                            for (auto r: repeat_contigs){
                                nc.canonical_repeats.back().push_back(r);
                            }
                            std::cout << "  nc.canonical_repeats 3 " <<  nc.canonical_repeats.size() << std::endl ;
                                      //
                            canonical.push_back(repeat_contigs);
                        } else {
                            outfile << "repeat not present " << std::endl;
                        }
                        outfile << std::endl;

                        local_repeat_contig_values.clear();
                    } else {
                        auto kci_node = kci.compute_compression_for_node(counter, 10, nc.index);
                        count += 1;
                        for (auto k:kci_node) {


                                outfile_csv << k << ", ";

                        }
                        nc.compressions[counter] =  kci_node[1]/kci_node[6];
                    }
                }
    }

outfile_csv << std::endl;
    std::cout << "nc.compressionssize " << nc.compressions.size() << " counter " << sg.nodes.size() << " count " << count  <<std::endl;
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
        std::cout << nc.compressions[r[0]] << " repeat: " << sg.oldnames[r[0]] << " "<<std::endl;
        for(auto c: r) {
            all_repeat_compressions.push_back(nc.compressions[c]);
            std::cout << nc.compressions[c] << " repeat: " << sg.oldnames[c] << " ";
            std::cout <<std::endl;

        }
    }

    std::cout << "calculated compression for " << nc.lib_name_r1 << " for " <<
              compressions.size() << " nodes, maps to more than 3 nodes of " << nc.canonical_repeats.size() << " repeats, of which " << resolved_repeat_count << "resolved. \n"<<
                     " including the in/out contigs,  " << all_repeat_contig_values.size() << " scored in total, of which  " << in_out_sane << " in/out compression sums are consistent and " << repeated_contig_sane << " repeated contig compressions are consistent with in/out suns"  << std::endl;


    return;

    //TODO:: do this properly!! make tedt datast with canibical repeats + read sets resolving!!!
        }



