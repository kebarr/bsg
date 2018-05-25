//
// Created by Katie Barr (EI) on 19/05/2018.
//

#include "UniqueContigIdentifier.h"

UniqueContigIdentifier::UniqueContigIdentifier(std::vector<std::string> filenames, uint64_t _max_mem): assembly_names(filenames), max_mem(_max_mem){

  /*for (int i=0; i < filenames.size() ; i++) {
      assembly_names[filenames[i]] = i;
  }

    max_mem= _max_mem;*/

    graphs.resize(assembly_names.size());
    contig_ids.resize(assembly_names.size());
    kmers_per_assembly.resize(assembly_names.size());
};

void UniqueContigIdentifier::LoadGraphs(){
    graphs.resize(assembly_names.size());
    for (int i =0;  i < assembly_names.size() ; i++){
        auto sg = BuildGraph(assembly_names[i]);
        graphs[i] = sg;
    }

};
bool UniqueContigIdentifier::UniqueKmersPerContig(int lib_index, sgNodeID_t node_id, double match_thresh=0.5) {
    return MatchingKmersPerContig(lib_index, node_id, unique_kmers_for_assemblies[lib_index], match_thresh);
} ;



bool UniqueContigIdentifier::CommonKmersPerContig(int lib_index, sgNodeID_t node_id, double match_thresh=0.5) {
    return MatchingKmersPerContig(lib_index, node_id, common_kmers_for_assemblies[lib_index],match_thresh);
} ;


bool UniqueContigIdentifier::MatchingKmersPerContig(int lib_index, sgNodeID_t node_id, std::vector<uint64_t > kmers, double match_thresh=0.5){
   auto node= graphs[lib_index].nodes[node_id];
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence, 31);
    nkmers.reserve(node.sequence.size());
    skf.create_kmers(nkmers);
    std::vector<uint64_t > found;
    int matches_required = nkmers.size()*match_thresh;
    int count = 0;
    for (auto kmer:kmers) {
        if (found.size() < nkmers.size() * 0.1 && count >= nkmers.size() * 0.5) {
            return false;//
        }
        if (std::find(nkmers.begin(), nkmers.end(), kmer) != nkmers.end()) {
            found.push_back(kmer);

            if (found.size() > matches_required){
                return  true;
            }
        }
    }

};

std::vector<sgNodeID_t > UniqueContigIdentifier::GetUniqueContigs(int lib_index){
    std::vector<sgNodeID_t > contigs;
    for (int n=0; graphs[lib_index].nodes.size(); n++){
        auto to_include = UniqueKmersPerContig(lib_index, n);
            if (to_include){
                contigs.push_back(n);
            }
    }
    std::cout << contigs.size() << " contigs from " << graphs[lib_index].filename << " with more than the thresghold level pf unique content "<< std::endl;
    return  contigs;

};


std::vector<std::vector<sgNodeID_t > > UniqueContigIdentifier::GetAllUniqueContigs(double kmer_thresh=0.05){
    GetAllUniqueKmers();
    std::vector<std::vector<sgNodeID_t > > res;
    for (int i = 0; i < graphs.size() ; i++){
        // if a certain proprtoon are unique...
        if (unique_kmers_for_assemblies[i].size() > kmers_per_assembly[i]*kmer_thresh){
                auto contigs = GetUniqueContigs(i);
res.push_back(contigs);
        } else {
            res.push_back({});
        }
    }

};


void UniqueContigIdentifier::GetAllUniqueKmers(){
    std::cout << " getting kmerd from "<< graphs.size() << " asssemblies " << std::endl;
    for (int i = 0; i < graphs.size() ; i++){
       auto k =  GetKmers(i);
        kmers_per_assembly[i] = k;
    }
    std::cout << " got" << all_kmers.size() <<  " kmerd from asssemblies " << std::endl;


    unique_kmers_for_assemblies.resize(graphs.size());
    int counter =0;
    int counter_common =0;
    for (auto k:all_kmers){
        if (k.second.size() == 1){
            unique_kmers_for_assemblies[*k.second.begin()].push_back(k.first);
            counter +=1;
        } else if (k.second.size() == graphs.size()){
            common_kmers_for_assemblies[*k.second.begin()].push_back(k.first);
            counter_common += 1;
        }
    }
    int mean = 0;
    int min = all_kmers.size();
    int max = 0;
    std::string min_a;
    std::string max_a;
    int count = 0;
    for (auto u:unique_kmers_for_assemblies){
        mean += u.size();
        if (u.size()> max){
            max = u.size();
            max_a = assembly_names[count];
        } else if (u.size() < min){
            min = u.size();
            min_a = assembly_names[count];
        }
        count+=1;
    }
    mean = mean/(double)unique_kmers_for_assemblies.size();
    std::cout << "found " << std::to_string(counter) << " kmers unique to one assembly, " << std::to_string(counter_common) << " ound in all, with average of  " << mean << " kmers unique to each assembly " << max << " max unique mkers, in " << max_a << " and min " << min << " unique mkers for "<< min_a <<std::endl;

    return unique_kmers_for_assemblies;
};


void UniqueContigIdentifier::WriteUniqueContentToFasta(int lib_index){};


SequenceGraph UniqueContigIdentifier::BuildGraph(std::string gfa_filename) {
    std::cout << " Building " << gfa_filename << std::endl;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::cout << " Built " << gfa_filename << " with " << sg.nodes.size() << std::endl;
    return sg;

};

int UniqueContigIdentifier::GetKmers(int i){
    auto sg = graphs[i];
    std::cout << " Adding kmers from  " << sg.filename << std::endl;
size_t  old_size = all_kmers.size();
    KmerCompressionIndex kci(sg, max_mem*1024L*1024L*1024L);
    kci.index_graph();
    auto gk = kci.graph_kmers;
int count = 0;
    for (auto k:gk){
       // if (k.count == 1){ - think in this sense its kmers unique to variety, not unique in graph
            count +=1;
            all_kmers[k.kmer].insert(i);
       // }
    }
    size_t  new_size = all_kmers.size() - old_size;

    std::cout << "Added " << std::to_string(new_size) << " novel kmers from " << gk.size() << " kmers in assembly" << std::endl;

    return  count;
};
