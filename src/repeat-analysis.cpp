//
// Created by Katie Barr (EI) on 06/04/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.h>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/GraphPartitioner.hpp>
#include "sglib/SequenceGraph.h"
#include "cxxopts.hpp"
#include <sglib/CompressionAnalyzer.h>
#include <sglib/RepeatAnalyzer.h>



int main(int argc, char * argv[]) {

    std::cout << "Welcome to kci-comparer" << std::endl << std::endl;
    sglib::OutputLogLevel = sglib::DEBUG;
    //std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    //std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;

    std::string gfa_filename, output_prefix, load_cidx, gfa_list, assembly_list, mode;
    std::vector<std::string> reads1, reads2, reads_type, dump_mapped, load_mapped, cidxreads1, cidxreads2, dump_cidx;
    bool stats_only = 0;
    uint64_t max_mem_gb = 4;
    int number_repeats = 10, in_min=3000, out_min=3000,  rep_min=3000;

    try {
        cxxopts::Options options("repeat-analyser", "Compare kmer compression indices from read and contig data");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file list", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Repeat analysis options")
                ("repeats", "number of repeats to find", cxxopts::value<int>(number_repeats))
                ("in_min", "minimum number of bases n in contigs",  cxxopts::value<int>(in_min))

                ("out_min", "minimum number of bases n in contigs",  cxxopts::value<int>(out_min))
                ("rep_min", "minimum number of bases n in contigs",  cxxopts::value<int>(rep_min));

        options.add_options("Compression Index Options")
                ("cidxread1", "compression index input reads, left",
                 cxxopts::value<std::vector<std::string>>(cidxreads1))
                ("cidxread2", "compression index input reads, right",
                 cxxopts::value<std::vector<std::string>>(cidxreads2))

                ("mode", "kci mode",
                 cxxopts::value<std::string>(mode))
                ("load_cidx", "load compression index filename", cxxopts::value<std::string>(load_cidx))
                ("dump_cidx", "dump compression index filename", cxxopts::value<std::vector<std::string>>(dump_cidx))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));


        auto result(options.parse(argc, argv));

        if (result.count("help")) {
            std::cout << options.help({"", "Compression Index Options"}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }


        if (result.count("cidxread1") != result.count("cidxread2")) {
            throw cxxopts::OptionException(" please specify cidxread1 and cidxread2 files in pairs");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::ofstream outfile;
    outfile.open(output_prefix, std::ofstream::out |std::ofstream::app);
    std::cout << "Executed command:" << std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;
    for (auto i = 0; i < argc; i++) std::cout << argv[i] << " ";
    std::cout << std::endl << std::endl;

    auto fasta_filename = gfa_filename.substr(0, gfa_filename.size() - 4) + ".fasta";
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    CompressionAnalyzer ca(sg, max_mem_gb, output_prefix +"_detailed");

    for (int lib = 0; lib < cidxreads1.size(); lib++) {
        if (dump_cidx.size() == 0) {
            ca.InitializeLib(cidxreads1[lib], cidxreads2[lib]);
        } else {
            ca.InitializeLib(cidxreads1[lib], cidxreads2[lib], dump_cidx[lib]);
        }

        outfile << "lib: " << lib << " " << cidxreads1[lib] << " " << cidxreads2[lib];

    }
}

