//
// Created by Katie Barr (EI) on 06/04/2018.
//

#include "RepeatAnalyzer.h"
#include <sglib/CompressionAnalyzer.h>


RepeatAnalyzer::RepeatAnalyzer(SequenceGraph &_sg, KmerCompressionIndex &_kci, std::string lib_name=""):sg(_sg), kci(_kci){

};

void RepeatAnalyzer::FindRepeats(std::string name_base="rep", int limit=-1, int rep_min = 3000, int in_min=1000,int out_min=1000){
    int lim= limit==-1? sg.nodes.size(): limit;

    for (sgNodeID_t counter = 1; counter < sg.nodes.size(); counter++) {
        if (repeats.size() < limit && sg.nodes[counter].sequence.size() < rep_min, sg.is_canonical_repeat(counter)){
            auto in = sg.get_bw_links(counter);

            std::vector<sgNodeID_t > in_ids = {};

            std::vector<sgNodeID_t > out_ids = {};

            for (auto c:in){
                if (sg.nodes[c.dest].sequence.size() < in_min){
                    continue;
                }
                in_ids.push_back(c.dest);
            }
            auto out = sg.get_fw_links(counter);
            for (auto c:out){
                if (sg.nodes[c.dest].sequence.size() < out_min){
                    continue;
                }
                out_ids.push_back(c.dest);
            }
            Repeat repeat{name_base, counter, in_ids, out_ids};
            this->repeats.push_back(repeat);
        }

    }
    read_sets_mapping_to_repeats.resize(repeats.size());
};


/*
 *
 * bj::
std::vector<std::unordered_set<uint64_t>> FlowFollower::get_distinctive_kmers(std::vector<sgNodeID_t> nodes) {
    std::unordered_set<uint64_t> seen_kmers,shared_kmers;
    std::vector<std::unordered_set<uint64_t>> distinctive_kmers;
    for (auto n:nodes){
        distinctive_kmers.emplace_back();
        StringKMerFactory skf(ws.sg.nodes[llabs(n)].sequence,31);
        std::vector<uint64_t> nkmers;
        nkmers.reserve(ws.sg.nodes[llabs(n)].sequence.size());
        skf.create_kmers(nkmers);
        for (auto x:nkmers) {
            if (seen_kmers.count(x) > 0) {
                shared_kmers.insert(x);
            } else {
                distinctive_kmers.back().insert(x);
                seen_kmers.insert(x);
            }
        }
    }
    for (auto &dk:distinctive_kmers) for (auto sk:shared_kmers) if (dk.count(sk)) dk.erase(sk);
    return distinctive_kmers;
}
*/


void RepeatAnalyzer::OutputRepeats(std::string fname, std::vector<size_t > to_include={}){

    std::ofstream outfile;
    outfile.open(fname);
    //int count =
    for (auto r:repeats) {
        outfile << ">" << r.name_base << "_" << r.repeated_contig << "\n"
                << sg.nodes[r.repeated_contig].sequence << "\n";
        for (auto in :r.in_contigs) {
        outfile << ">" << r.name_base << "_" << in << "\n"
                << sg.nodes[in].sequence << "\n";
    }
        for (auto out :r.in_contigs) {
            outfile << ">" << r.name_base << "_" << out << "\n"
                    << sg.nodes[out].sequence << "\n";
        }

    }

};

std::vector<double> RepeatAnalyzer::compressions_for_read_set(double (*compression_function)(sgNodeID_t, KmerCompressionIndex &) ){
    int index=0;
    int repeats_mapped=0;
    std::vector<double> compressions;
    std::cout << "Computing compressions for " << kci.current_lib << std::endl;
    for (auto r:repeats){
        RepeatReduction(r);
        RepeatCompressions rc("kci.compute_kcov_for_node", compression_function,  r.repeated_contig,kci, r.in_contigs, r.out_contigs);
        r.rc["kci.compute_kcov_for_node"].push_back(rc);
        if (rc.repeat_contained && rc.in_out_agree && rc.repeat_coverage_sane) {read_sets_mapping_to_repeats[index].push_back(kci.current_lib);
        repeats_mapped+=1;
            for (auto c:rc.in_compressions){
                compressions.push_back(c);

            }
            compressions.push_back(rc.repeat_compression);
            for (auto c:rc.out_compressions){
                compressions.push_back(c);

            }
        };


    }
    std::cout << "Read set " << kci.current_lib  << " mapped to " << repeats_mapped << " of " << repeats.size() << " repeats in  " << sg.filename<< std::endl;

    return  compressions;
};

std::vector<std::vector<uint64_t> > RepeatAnalyzer::ConvertParallelContigsDistinctKmers(std::vector<sgNodeID_t > contigs){
    std::vector<std::vector<uint64_t> > node_kmers;
    size_t to_compare = sg.nodes[contigs[0]].sequence.size();
    for (int c = 0; c< contigs.size(); c++){
        if (sg.nodes[contigs[c]].sequence.size() < to_compare){
            to_compare = sg.nodes[contigs[c]].sequence.size();
        }
    }
    for (auto c:contigs) {
        StringKMerFactory skf(sg.nodes[llabs(c)].sequence,31);
        std::vector<uint64_t> nkmers;
        node_kmers.push_back(nkmers);
        nkmers.reserve(sg.nodes[llabs(c)].sequence.size());
        skf.create_kmers(nkmers);
        node_kmers.push_back(nkmers);
    }
    return node_kmers;
};

// thiRepeatAnalyzer::s can't return an actual repeat as reduced contigs not in graph, so represent result on repeat struct, each in/out contig shoulxd just be a sequence of distinct kmers + positions in contig
void RepeatAnalyzer::RepeatReduction(Repeat repeat){
 /*   std::vector<sgNodeID_t > nodes_in;
    std::vector<sgNodeID_t > nodes_out;
    std::vector<std::vector<uint64_t> > node_kmers;
    auto in_kmers = ConvertParallelContigsDistinctKmers(repeat.in_contigs);
    auto out_kmers = ConvertParallelContigsDistinctKmers(repeat.out_contigs);

    for (auto in:repeat.in_contigs) {nodes_in.push_back(in);
        std::vector<uint64_t> nkmers;
        StringKMerFactory skf(sg.nodes[llabs(in)].sequence,31);
        node_kmers.push_back(nkmers);
        nkmers.reserve(sg.nodes[llabs(in)].sequence.size());
        skf.create_kmers(nkmers);}
    for (auto out:repeat.out_contigs) {nodes_out.push_back(out);}

    //startig from closrst to repeat contig, include only kmers that differ
    std::unordered_set<uint64_t> seen_kmers,shared_kmers;
    size_t to_cpomare_in = sg.nodes[nodes_in[0]].sequence.size();
    for (int c = 0; c< nodes_in.size(); c++){
        if (sg.nodes[nodes_in[c]].sequence.size() < to_cpomare_in){
            to_cpomare_in = sg.nodes[nodes_in[c]].sequence.size();
        }
    }
    std::vector<std::vector<uint64_t> > node_kmers;
    std::vector<std::unordered_set<uint64_t>> distinctive_kmers;
    for (auto n:nodes){
        distinctive_kmers.emplace_back();
        StringKMerFactory skf(sg.nodes[llabs(n)].sequence,31);
        std::vector<uint64_t> nkmers;
        node_kmers.push_back(nkmers);
        nkmers.reserve(sg.nodes[llabs(n)].sequence.size());
        skf.create_kmers(nkmers);
        for (auto x:nkmers) {
            if (seen_kmers.count(x) > 0) {
                shared_kmers.insert(x);
            } else {
                distinctive_kmers.back().insert(x);
                seen_kmers.insert(x);
            }
        }
    }
    std::vector<std::vector<uint64_t > > in_contigs_distinct_kmers;
    in_contigs_distinct_kmers.resize(repeat.in_contigs.size());
    for (int i=0; i < to_cpomare_in; i++){

    } */
};