//
// Created by Katie Barr (EI) on 19/05/2018.
//

#include "UniqueContigIdentifier.h"

UniqueContigIdentifier::UniqueContigIdentifier(std::vector<std::string> filenames, uint64_t _max_mem): assembly_names(filenames), max_mem(_max_mem){

  /*for (int i=0; i < filenames.size() ; i++) {
      assembly_names[filenames[i]] = i;
  }

    max_mem= _max_mem;*/

    LoadGraphs();
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
    return MatchingKmersPerContig(lib_index, node_id, common_kmers_for_assemblies,match_thresh);
} ;

void UniqueContigIdentifier::FindCoreGenome(int base_lib){
    auto uniq = GetUniqueContigs(base_lib);
    auto next = uniq[0];
    int ind=1;
    std::vector<sgNodeID_t > common1;
    std::vector<sgNodeID_t > common2;
    for (sgNodeID_t  n=0 ; n < graphs[base_lib].nodes.size() ; n++){
        if (n != next){
            common1.push_back(n);
        } else {
            next = uniq[ind];
            ind +=1;
        }
        if (CommonKmersPerContig(base_lib, n)){
            common2.push_back(n);
        }
    }
    std::cout << common1.size() << " core contigs with method 1 " << common2.size() << " core contigs with method 2 " << std::endl;


}


bool UniqueContigIdentifier::MatchingKmersPerContig(int lib_index, sgNodeID_t node_id, std::vector<uint64_t > kmers, double match_thresh=0.5){
   auto node= graphs[lib_index].nodes[node_id];
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence, 31);
    nkmers.reserve(node.sequence.size());
    skf.create_kmers(nkmers);
    std::vector<uint64_t > found;
    auto matches_required = nkmers.size()*match_thresh;
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


void UniqueContigIdentifier::GetAllUniqueContigs(double kmer_thresh=0.05){
    GetAllUniqueKmers();
    std::vector<std::vector<sgNodeID_t > > res;
    for (int i = 0; i < graphs.size() ; i++){
        // if a certain proprtoon are unique...
        if (unique_kmers_for_assemblies[i].size() > kmers_per_assembly[i]*kmer_thresh){
                auto contigs = GetUniqueContigs(i);
res.push_back(contigs);
            std::cout << contigs.size() << "  unique contigs for " << assembly_names[i] << " which cntains " << kmers_per_assembly[i] << "kmers " << graphs[i].nodes.size() << " contigs " << unique_kmers_for_assemblies[i].size() << " unique kmers " <<  std::endl;
        } else {
            res.emplace_back({});
        }
    }
    this->unique_contigs = res;

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
    for (auto k:all_kmers){// kmer: assemblies included in
        if (k.second.size() == 1){
            unique_kmers_for_assemblies[*k.second.begin()].push_back(k.first);
            counter +=1;
        } else if (k.second.size() == graphs.size() || k.second.size() == graphs.size() -1){
            for (auto s:k.second) {
                common_kmers_for_assemblies.push_back(k.first);
            }
            counter_common += 1;
        }
    }
    double mean = 0;
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

    std::cout << std::to_string(common_kmers_for_assemblies.size()) << " kmers found in almost all assemblies " << std::endl;
};


void UniqueContigIdentifier::WriteUniqueContentToFasta(int lib_index, std::string fname){
    std::ofstream out;
    out.open(fname);
    for (auto c:unique_contigs[lib_index]){
        out << ">" <<  graphs[lib_index].oldnames[c]<< std::endl << graphs[lib_index].nodes[c].sequence<< std::endl;

    }
};


SequenceGraph UniqueContigIdentifier::BuildGraph(std::string gfa_filename) {
    std::cout << " Building " << gfa_filename << std::endl;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::cout << " Built " << gfa_filename << " with " << sg.nodes.size() << std::endl;
    return sg;

};

size_t UniqueContigIdentifier::GetKmers(int i){
    auto sg = graphs[i];
    std::cout << " Adding kmers from  " << sg.filename << std::endl;
size_t  old_size = all_kmers.size();
    KmerCompressionIndex kci(sg, max_mem*1024L*1024L*1024L);
    kci.index_graph();
    auto gk = kci.graph_kmers;
    for (auto k:gk){
       // if (k.count == 1){ - think in this sense its kmers unique to variety, not unique in graph
            all_kmers[k.kmer].insert(i);
       // }
    }
    size_t  new_size = all_kmers.size() - old_size;

    std::cout << "Added " << std::to_string(new_size) << " novel kmers from " << gk.size() << " kmers in assembly" << std::endl;

    return  gk.size();
};
