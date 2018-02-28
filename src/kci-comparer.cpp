#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.h>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/GraphPartitioner.hpp>
#include "sglib/SequenceGraph.h"
#include "cxxopts.hpp"


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

    sglib::OutputLogLevel = sglib::DEBUG;
    std::cout << "Welcome to kci-comparer"<<std::endl<<std::endl;
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
            kci.start_new_count();
            kci.add_counts_from_file(cidxreads1[lib]);
            kci.add_counts_from_file(cidxreads2[lib]);
            kci.compute_compression_stats(lib);
            kci.dump_histogram(output_prefix + "_" + std::to_string(lib) + ".csv", lib);
            std::cout << "Counted reads for lib " << lib << " \n";
            kci_assembly2 << "lib: " << lib << " " << cidxreads1[lib] << " " << cidxreads2[lib];
            for (sgNodeID_t counter = 0; counter < sg.nodes.size(); counter++) {

                if (sg.is_canonical_repeat(counter)) {
                    auto bw = sg.get_bw_links(counter);
                    auto fw = sg.get_fw_links(counter);
                    // percent present/absent doesn't do it - or not obviously
                    for (auto b: bw){
                        auto kci_node = kci.compute_compression_for_node(b.dest, 10, lib);
                        kci_assembly << kci_node << ", ";
                        auto ind = b.dest > 0 ? b.dest : -b.dest;
                        kci_assembly2 << sg.oldnames[ind] << ": " << kci_node << ", " << sg.nodes[ind].sequence.size() << ", ";
                    }
                    auto kci_node = kci.compute_compression_for_node(counter, 10, lib);
                    kci_assembly << kci_node << ", ";
                    kci_assembly2 << sg.oldnames[counter] << ": " << kci_node << ", "<< sg.nodes[counter].sequence.size() << ", ";
                    for (auto f: fw){
                        auto kci_node = kci.compute_compression_for_node(f.dest, 10, lib);
                        kci_assembly << kci_node << ", ";
                        auto ind = f.dest > 0? f.dest : -f.dest;
                        kci_assembly2 << sg.oldnames[ind] << ": " << kci_node << ", " << sg.nodes[ind].sequence.size() << ", ";
                    }
                    kci_assembly2 << std::endl;
                    count += 1;

                } else {
                    kci_assembly << " ,";
                }

            }
            kci_assembly << std::endl;

            std::cout << "calculated compression for " << lib << " for " << count << std::endl;
        }
    }

}
