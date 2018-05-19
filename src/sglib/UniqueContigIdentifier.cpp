//
// Created by Katie Barr (EI) on 19/05/2018.
//

#include "UniqueContigIdentifier.h"

UniqueContigIdentifier::UniqueContigIdentifier(std::vector<std::string> filenames, uint64_t _max_mem){

  for (int i=0; i < filenames.size() ; i++) {
      assembly_names[filenames[i]] = i;
  }

    max_mem= _max_mem;

};

void UniqueContigIdentifier::LoadGraphs(){
    graphs.resize(assembly_names.size());
    for (auto f:assembly_names){
        auto sg = BuildGraph(f.first);
        graphs[f.second] = sg;
    }

};

std::vector<sgNodeID_t > UniqueContigIdentifier::GetUniqueContigs(int lib_index){}
;
std::vector<sgNodeID_t > UniqueContigIdentifier::GetAllUniqueContigs(){

};


std::vector<sgNodeID_t > UniqueContigIdentifier::GetAllUniqueKmers(){
    std::cout << " getting kmerd from asssemblies " << std::endl;
    for (auto g:graphs){
        GetKmers(g);
    }
    std::cout << " got" << all_kmers.size() <<  " kmerd from asssemblies " << std::endl;

    std::vector<std::vector< uint64_t > > unique_kmers_for_assemblies;
    unique_kmers_for_assemblies.resize(graphs.size());
    int counter =0;
    for (auto k:all_kmers){
        if (k.second.size() == 1){
            unique_kmers_for_assemblies[assembly_names[*k.second.begin()]].push_back(k.first);
            counter +=1;
        }
    }
    int mean = 0;
    int min = all_kmers.size();
    int max = 0;
    int min_a;
    int max_a;
    int count = 0;
    for (auto u:unique_kmers_for_assemblies){
        mean += u.size();
        if (u.size()> max){
            max = u.size();
            max_a = count;
        } else if (u.size() < min){
            min = u.size();
            min_a = count;
        }
        count+=1;
    }
    mean = mean/(double)unique_kmers_for_assemblies.size();
    std::cout << "found " << std::to_string(counter) << " kmers unique to one assembly, with average of  " << mean << " kmers unique to each assembly " <<  <<std::endl;
};


void UniqueContigIdentifier::WriteUniqueContentToFasta(int lib_index){};


SequenceGraph UniqueContigIdentifier::BuildGraph(std::string gfa_filename) {
    std::cout << " Building " << gfa_filename << std::endl;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::cout << " Built " << gfa_filename << " with " << sg.nodes.size() << std::endl;
    return; SequenceGraph;

};

void UniqueContigIdentifier::GetKmers(SequenceGraph & sg){
    std::cout << " Adding kmers from  " << sg.filename << std::endl;
size_t  old_size = all_kmers.size();
    KmerCompressionIndex kci(sg, max_mem*1024L*1024L*1024L)
    kci.index_graph();
    auto gk = kci.graph_kmers;
    std::vector<uint64_t > kmers;
    for (auto k:gk){
       // if (k.count == 1){ - think in this sense its kmers unique to variety, not unique in graph
            kmers.push_back(k.kmer);
            all_kmers[k.kmer].insert(sg.filename);
       // }
    }
    size_t  new_size = all_kmers.size() - old_size;

    std::cout << "Added " << std::to_string(new_size) << " novel kmers from " << gk.size() << " kmers in assembly" << std::endl;
};
