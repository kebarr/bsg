#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.h>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/GraphPartitioner.hpp>
#include "sglib/SequenceGraph.h"
#include "cxxopts.hpp"



std::vector<double> analyse_repeat(std::vector<double> repeat_compressions, double tolerance=0.95, double diff_threshold=10){

    // see if repeat phased by reads- if compressions on each side pair to be reasonably close
    // sum of compressions on each side should be a multiple of middle compression

    auto in_sum = repeat_compressions[0] + repeat_compressions[1];
    auto out_sum = repeat_compressions[3] + repeat_compressions[4];
    double in_out_sane = 0;
    double repeat_count_sane = 0;

    std::vector<double> res;
    if (tolerance*in_sum < out_sum < (1-tolerance)*in_sum){
        in_out_sane = 1;
        // not sure this actually works for way i\m calculating 'compression'
        // contig repeated 5 timea should have 5*kmers in reads than average, and 5*reads going in, split in a sane way- i/e. shouldn't be 1 kmer on one in, 100 on other in, then 50/50 out
        if (std::abs(in_sum - repeat_compressions[2]) < tolerance && std::abs(out_sum - repeat_compressions[2]) < tolerance){
            repeat_count_sane = 1;
        }
        // in this ase check if resolves repeat, find out closest to in and see if close enough to call
        auto pairs = repeat_compressions[3]-repeat_compressions[0] < repeat_compressions[4]-repeat_compressions[0] ? std::make_pair(3, 4) : std::make_pair(4, 3);
        if (std::abs(repeat_compressions[0] - repeat_compressions[std::get<0>(pairs)]) < std::abs((repeat_compressions[0] - repeat_compressions[std::get<1>(pairs)]*diff_threshold)) && std::abs(repeat_compressions[1] - repeat_compressions[std::get<1>(pairs)]) < std::abs((repeat_compressions[0] - repeat_compressions[std::get<1>(pairs)]*diff_threshold))){
            // then accprding to this arbitrary heiristic, we resolve to get 0 with pair 0
            res.push_back(repeat_compressions[std::get<0>(pairs)]);
            res.push_back(repeat_compressions[0] + repeat_compressions[std::get<0>(pairs)]);

            res.push_back(repeat_compressions[std::get<1>(pairs)]);
            res.push_back(repeat_compressions[1] + repeat_compressions[std::get<1>(pairs)]);

        } else {
            // if its not resolved its more useful to know how different they were
            res.push_back(0);
            res.push_back(repeat_compressions[0] - repeat_compressions[std::get<0>(pairs)]);
            res.push_back(0);
            res.push_back(repeat_compressions[1] - repeat_compressions[std::get<1>(pairs)]);

        }
        res.push_back(in_out_sane);
        res.push_back(repeat_count_sane);
    }
    return  res;

}


void output_kci_for_assembly(std::string gfa_name,std::string assembly_name, std::string cidxread1, std::string cidxread2, uint64_t max_mem_gb, std::string dump_cidx=""){
    auto fasta_filename=gfa_name.substr(0,gfa_name.size()-4)+".fasta";
    SequenceGraph sg;
    sg.load_from_gfa(gfa_name);

std::cout<<std::endl<<"=== Loading reads compression index ==="<<std::endl;
//compression index
KmerCompressionIndex kci(sg,max_mem_gb*1024L*1024L*1024L);

    kci.index_graph();
        kci.start_new_count();
        kci.add_counts_from_file(cidxread1);
        kci.add_counts_from_file(cidxread2);

if (dump_cidx!=""){
    kci.save_to_disk(dump_cidx);
}

if (kci.read_counts.size()>0) {
    kci.compute_compression_stats();
    kci.dump_histogram(assembly_name + "kci_histogram.csv");
}
std::string outname = assembly_name+ "_kci_per_node.csv";
std::ofstream kci_assembly(outname);
for (size_t counter = 0; counter < sg.nodes.size(); counter++){
    kci_assembly<< sg.oldnames[counter] << ", ";
}
kci_assembly<< std::endl;
for (sgNodeID_t counter = 0; counter < sg.nodes.size(); counter++){
    auto kci_node = kci.compute_compression_for_node(counter);
    kci_assembly << kci_node <<", ";

}
kci_assembly << std::endl;
}

void output_kci_for_read_set(SequenceGraph &sg, KmerCompressionIndex & kci, std::string assembly_name, std::string cidxread1, std::string cidxread2,std::ofstream kci_assembly, std::string dump_cidx=""){

    kci.start_new_count();
    kci.add_counts_from_file(cidxread1);
    kci.add_counts_from_file(cidxread2);

    if (dump_cidx!=""){
        kci.save_to_disk(dump_cidx);
    }

    if (kci.read_counts.size()>0) {
        kci.compute_compression_stats();
        kci.dump_histogram(assembly_name + "kci_histogram.csv");
    }
    for (sgNodeID_t counter = 0; counter < sg.nodes.size(); counter++){
        auto kci_node = kci.compute_compression_for_node(counter);
        kci_assembly << kci_node <<", ";

    }
    kci_assembly << std::endl;
}


int main(int argc, char * argv[]) {

    std::cout << "Welcome to kci-comparer"<<std::endl<<std::endl;
    sglib::OutputLogLevel = sglib::DEBUG;
    //std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    //std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;

    std::string gfa_filename,output_prefix, load_cidx, dump_cidx, gfa_list, assembly_list;
    std::vector<std::string> reads1,reads2,reads_type, dump_mapped, load_mapped,cidxreads1,cidxreads2;
    bool stats_only=0;
    uint64_t max_mem_gb=4;

    try
    {
        cxxopts::Options options("kci-comparer", "Compare kmer compression indices from read and contig data");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file list", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Compression Index Options")
                ("cidxread1", "compression index input reads, left", cxxopts::value<std::vector<std::string>>(cidxreads1))
                ("cidxread2", "compression index input reads, right", cxxopts::value<std::vector<std::string>>(cidxreads2))
                ("load_cidx", "load compression index filename", cxxopts::value<std::string>(load_cidx))
                ("dump_cidx", "dump compression index filename", cxxopts::value<std::string>(dump_cidx))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb))
                ("gfa_list", "list of gfas for compairson", cxxopts::value<std::string>(gfa_list));


        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Compression Index Options"}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }


        if ( result.count("cidxread1")!=result.count("cidxread2")){
            throw cxxopts::OptionException(" please specify cidxread1 and cidxread2 files in pairs");
        }




    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;
    if (gfa_list != "") {
        std::ifstream infile(gfa_list);
        std::string line;
        std::string fields[2];
        while (std::getline(infile, line)) {
            std::istringstream(line) >> fields[0] >> fields[1];
            std::cout << "calculating compression for: " << fields[0] << " from " << fields[1] << std::endl;
            output_kci_for_assembly(fields[0], fields[1], cidxreads1[0], cidxreads2[0], max_mem_gb, dump_cidx);

        }
    } else {

        auto fasta_filename = gfa_filename.substr(0, gfa_filename.size() - 4) + ".fasta";
        SequenceGraph sg;
        sg.load_from_gfa(gfa_filename);

        std::cout << std::endl << "=== Loading reads compression index ===" << std::endl;
//compression index
        KmerCompressionIndex kci(sg, max_mem_gb * 1024L * 1024L * 1024L);

        kci.index_graph();



        std::ofstream kci_assembly(output_prefix + "_kcis.csv");


        std::ofstream kci_assembly2(output_prefix + "_kcis_repeats.csv");
        for (size_t counter = 0; counter < sg.nodes.size(); counter++) {
            kci_assembly << sg.oldnames[counter] << ", ";
        }
        kci_assembly << std::endl;

        for(int lib=0;lib<cidxreads1.size();lib++) {
            int count = 0;
            int mapped_repeat_count = 0;
            int resolved_repeat_count = 0;

            int sum = 0;
            double in_out_sane = 0
            double repeated_contig_sane = 0;
            std::vector<double> compressions;
            std::vector<double> repeat_vals;
            std::vector<double> repeat_contig_values;
            kci.start_new_count();
            kci.add_counts_from_file(cidxreads1[lib]);
            kci.add_counts_from_file(cidxreads2[lib]);
            kci.compute_compression_stats();
            kci.dump_histogram(output_prefix + "_" + std::to_string(lib) + ".csv");
            std::cout << "Counted reads for lib " << lib << " \n";
            kci_assembly2 << "lib: " << lib << " " << cidxreads1[lib] << " " << cidxreads2[lib];
            for (sgNodeID_t counter = 0; counter < sg.nodes.size(); counter++) {

                if (sg.is_canonical_repeat(counter)) {
                    int nonzeros = 0;
                    auto bw = sg.get_bw_links(counter);
                    auto fw = sg.get_fw_links(counter);
                    // percent present/absent doesn't do it - or not obviously, now dropped unique kmer requirement
                    for (auto b: bw){
                        auto kci_node = kci.compute_compression_for_node(b.dest, 10, lib);
                        kci_assembly << kci_node << ", ";
                        auto ind = b.dest > 0 ? b.dest : -b.dest;
                        kci_assembly2 << sg.oldnames[ind] << ": " << kci_node << ", " << sg.nodes[ind].sequence.size() << ", ";
                        sum +=kci_node;
                        compressions.push_back(kci_node);
                        repeat_vals.push_back(kci_node);
                        if (kci_node > 0 ) nonzeros +=1;
                        kci_node = 0;

                    }
                    auto kci_node = kci.compute_compression_for_node(counter, 10, lib);
                    kci_assembly << kci_node << ", ";
                    sum +=kci_node;
                    compressions.push_back(kci_node);
                    repeat_contig_values.push_back(kci_node);
                    repeat_vals.push_back(kci_node);
                    if (kci_node > 0 ) nonzeros +=1;
                    kci_node = 0;


                    kci_assembly2 << sg.oldnames[counter] << ": " << kci_node << ", "<< sg.nodes[counter].sequence.size() << ", ";
                    for (auto f: fw){
                        auto kci_node = kci.compute_compression_for_node(f.dest, 10, lib);
                        kci_assembly << kci_node << ", ";
                        auto ind = f.dest > 0? f.dest : -f.dest;
                        kci_assembly2 << sg.oldnames[ind] << ": " << kci_node << ", " << sg.nodes[ind].sequence.size() << ", ";
                        sum +=kci_node;
                        compressions.push_back(kci_node);
                        repeat_vals.push_back(kci_node);
                        if (kci_node > 0 ) nonzeros +=1;

                        kci_node = 0;


                    }
                    kci_assembly2 << std::endl;
                    count += 1;
                    if (nonzeros >=3 ) {
                        mapped_repeat_count += 1;
                        auto res = analyse_repeat(repeat_vals);
                        in_out_sane += res[4];
                        repeated_contig_sane += res[5];

                        if (res[0] != 0) {
                            kci_assembly2 <<  "Resolved: ";
                            resolved_repeat_count += 1;
                        }

                        for (auto r:res ){
                            kci_assembly2 << r << ", ";
                        }
                    } else {
                        kci_assembly2 << "repeat not present " << std::endl;
                    }
                    kci_assembly2 << std::endl;

                    repeat_vals.clear();
                } else {
                    kci_assembly << " ,";
                }

            }
            kci_assembly << std::endl;
            double sd;
            double mean = sum/compressions.size();
            for (auto c:compressions){
                sd += (c-mean)*(c-mean);
            }
            sd = (std::pow(sd, 0.5))/(compressions.size()-1);
            auto max = std::max_element(compressions.begin(), compressions.end());
            auto min = std::min_element(compressions.begin(), compressions.end());

            std::cout << "calculated compression for " << lib << " for " <<
                                                                         compressions.size() << " nodes, maps to " << mapped_repeat_count << " repeats, of which " << resolved_repeat_count << "resolved. \n"
                    " calxulated scores for  " << repeat_contig_values.size() << " repeated contigs, including the in/out contigs,  " << repeat_vals.size() << " scored in total, of which  " << in_out_sane << " in/out compression sums are consistent and " << repeated_contig_sane << " repeated contig compressions are consistent with in/out suns"  << std::endl;

            std::cout << "stats for all contigs in cononical repeats, min  " << min* << " max "
            double sd2;
            double mean2 = sum/repeat_contig_values.size();
            for (auto c:repeat_contig_values){
                sd2 += (c-mean2)*(c-mean2);
            }
            sd2 = (std::pow(sd2, 0.5))/(repeat_contig_values.size()-1);
            auto max2 = std::max_element(repeat_contig_values.begin(), repeat_contig_values.end());
            auto min2 = std::min_element(repeat_contig_values.begin(), repeat_contig_values.end());
        }
    }

}
