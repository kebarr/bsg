//
// Created by Katie Barr (EI) on 06/04/2018.
//

#include "RepeatAnalyzer.h"
#include <sglib/CompressionAnalyzer.h>


RepeatAnalyzer::RepeatAnalyzer(SequenceGraph &_sg, std::string lib_name=""):sg(_sg){

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

void RepeatAnalyzer::RepeatReduction(std::vector<sgNodeID_t > nodes){
    //startig from closrst to repeat contig, include only kmers that differ
    std::unordered_set<uint64_t> seen_kmers,shared_kmers;
    size_t to_cpomare = sg.nodes[nodes[0]].sequence.size();
    for (int c = 1; c< nodes.size(); c++){
        if (sg.nodes[c].sequence.size() < to_cpomare){
            to_cpomare = sg.nodes[c].sequence.size();
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
    for (int i=0; i < to_cpomare; i++){

    }
};