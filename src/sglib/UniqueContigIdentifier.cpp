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
    variants.resize(assembly_names.size());
    for (int i =0;  i < assembly_names.size() ; i++){
        auto sg = BuildGraph(assembly_names[i]);
        Variant var;
        var.sg = sg;
        variants[i] = var;
    }

};
bool UniqueContigIdentifier::UniqueKmersPerContig(int lib_index, sgNodeID_t node_id, double match_thresh=0.5) {
    return MatchingKmersPerContig(lib_index, node_id, unique_kmers_for_assemblies[lib_index], match_thresh);
} ;



bool UniqueContigIdentifier::CommonKmersPerContig(int lib_index, sgNodeID_t node_id, double match_thresh=0.5) {
    return MatchingKmersPerContig(lib_index, node_id, common_kmers_for_assemblies,match_thresh);
} ;


bool UniqueContigIdentifier::NonCoreKmersPerContig(int lib_index, sgNodeID_t node_id, double match_thresh=0.5) {
    return MatchingKmersPerContig(lib_index, node_id, non_core_kmers_for_assemblies[lib_index],match_thresh);
} ;

void UniqueContigIdentifier::FindCoreGenome(int base_lib){
    GetUniqueContigs(base_lib);
    auto uniq = variants[base_lib].unique_contigs;
    auto non_core = variants[base_lib].non_core_contigs;
    auto next = uniq[0];
    auto nex2 = non_core[0];

    int ind=1;
    std::vector<sgNodeID_t > common1;
    std::vector<sgNodeID_t > common3;

    std::vector<sgNodeID_t > common2;
    for (sgNodeID_t  n=0 ; n < variants[base_lib].sg.nodes.size() ; n++){
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
   auto node= variants[lib_index].sg.nodes[node_id];
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

void UniqueContigIdentifier::SummarizeVars(){
    for (int i = 0; i < variants.size(); i++){
        auto v1 = variants[i];
        std::vector<uint64_t > no_core_v1;
        std::set_difference(v1.kmer_total.begin(), v1.kmer_total.end(), v1.common_kmers.begin(), v1.common_kmers.end(),std::back_inserter(no_core_v1));

        for (int j = i +1 ; j < variants.size(); j++){
            auto v2 = variants[j];
            std::vector<uint64_t > no_core_v2;

            std::set_difference(v2.kmer_total.begin(), v2.kmer_total.end(), v2.common_kmers.begin(), v2.common_kmers.end(), std::back_inserter(no_core_v2));

            std::vector<uint64_t > shared_kmers_all;
            std::set_intersection(v1.kmer_total.begin(), v1.kmer_total.end(), v2.kmer_total.begin(), v2.kmer_total.end(), std::back_inserter(shared_kmers_all));

            std::vector<uint64_t > shared_kmers_non_core;
            std::set_intersection(no_core_v1.begin(), no_core_v1.end(), no_core_v2.begin(), no_core_v2.end(), std::back_inserter(shared_kmers_non_core));

            v1.shared_variation_all[j] = shared_kmers_all.size()/(double) v1.kmer_total.size();
            v2.shared_variation_all[i] = shared_kmers_all.size()/(double) v2.kmer_total.size();

            v1.shared_variation_non_core[j] = shared_kmers_non_core.size()/(double) v1.kmer_total.size();
            v2.shared_variation_non_core[i] = shared_kmers_non_core.size()/(double) v2.kmer_total.size();
            this->variants[i] = v1;
            this-> variants[j] = v2;
        }
        std::cout << " Var: " << v1.sg.filename << " kmers: " << v1.kmer_total.size() << " unique: " << v1.unique_kmers.size() << " common: " << v1.common_kmers.size() << " non core: " << v1.non_core_kmers.size() << " ";
     int l=0;
        for (auto c: v1.shared_variation_non_core) {
            std::cout<< l << ": " << c << " ";
            l++;
        }
        int k=0;
        for (auto c: v1.shared_variation_all) {
            std::cout<< k << ": " << c << " ";
            k++;
        }

    }
}

void UniqueContigIdentifier::GetUniqueContigs(int lib_index){
    std::vector<sgNodeID_t > contigs;
    for (int n=0; variants[lib_index].sg.nodes.size(); n++){

                auto common = CommonKmersPerContig(lib_index, n);
                if (common) {
                    variants[lib_index].common_contigs.insert(n);
                } else {
                    variants[lib_index].non_core_contigs.push_back(n);
                    auto to_include = UniqueKmersPerContig(lib_index, n);
                    if (to_include){
                        contigs.push_back(n);
                        variants[lib_index].unique_contigs.push_back(n);
                    }
                }


    }
    std::cout << contigs.size() << " contigs from " << variants[lib_index].sg.filename << " with more than the thresghold level pf unique content "<< std::endl;

};


void UniqueContigIdentifier::GetAllUniqueContigs(double kmer_thresh=0.05){
    GetAllUniqueKmers();
    std::vector<std::vector<sgNodeID_t > > res;
    for (int i = 0; i < variants.size() ; i++){
        // if a certain proprtoon are unique...
        if (unique_kmers_for_assemblies[i].size() > kmers_per_assembly[i]*kmer_thresh){
                GetUniqueContigs(i);
res.push_back(variants[i].unique_contigs);
            std::cout << variants[i].unique_contigs.size() << "  unique contigs for " << assembly_names[i] << " which cntains " << kmers_per_assembly[i] << "kmers " << variants[i].sg.nodes.size() << " contigs " << unique_kmers_for_assemblies[i].size() << " unique kmers " <<  std::endl;
        } else {
            res.emplace_back();
        }
    }
    this->unique_contigs = res;

};


void UniqueContigIdentifier::GetAllUniqueKmers(){
    std::cout << " getting kmerd from "<< variants.size() << " asssemblies " << std::endl;
    for (int i = 0; i < variants.size() ; i++){
       auto k =  GetKmers(i);
        kmers_per_assembly[i] = k;
    }
    std::cout << " got" << all_kmers.size() <<  " kmerd from asssemblies " << std::endl;


    unique_kmers_for_assemblies.resize(variants.size());
    int counter =0;
    int counter_common =0;
    for (auto k:all_kmers){// kmer: assemblies included in
        if (k.second.size() == 1){
            unique_kmers_for_assemblies[*k.second.begin()].push_back(k.first);
            variants[*k.second.begin()].unique_kmers.insert(k.first);
            counter +=1;
        } else if (k.second.size() == variants.size() || k.second.size() == variants.size() -1){
            for (auto s:k.second) {
                variants[s].common_kmers.insert(k.first);
            }
            common_kmers_for_assemblies.push_back(k.first);

            counter_common += 1;
        }else if (k.second.size() > 0){
            for (auto s:k.second) {
                variants[s].non_core_kmers.insert(k.first);

                non_core_kmers_for_assemblies[s].push_back(k.first);
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

void UniqueContigIdentifier::WriteAllUniqueContentToFasta(){
    for (int i = 0; i < variants.size(); i++){
        auto filename = variants[i].sg.filename.substr(0, variants[i].sg.filename.size() - 4) + "unique.fasta";

        WriteUniqueContentToFasta(i, filename);
    }
}


void UniqueContigIdentifier::WriteAllNonCoreContentToFasta(){
    for (int i = 0; i < variants.size(); i++){
        auto filename = variants[i].sg.filename.substr(0, variants[i].sg.filename.size() - 4) + "unique.fasta";

        WriteNonCoreContentToFasta(i, filename);
    }
}

void UniqueContigIdentifier::WriteNonCoreContentToFasta(int lib_index, std::string fname){
    std::ofstream out;
    out.open(fname);
    for (auto c:variants[lib_index].non_core_contigs){
        out << ">" <<  variants[lib_index].sg.oldnames[c]<< std::endl << variants[lib_index].sg.nodes[c].sequence<< std::endl;

    }
};


void UniqueContigIdentifier::WriteUniqueContentToFasta(int lib_index, std::string fname){
    std::ofstream out;
    out.open(fname);
    for (auto c:variants[lib_index].unique_contigs){
        out << ">" <<  variants[lib_index].sg.oldnames[c]<< std::endl << variants[lib_index].sg.nodes[c].sequence<< std::endl;

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
    auto sg = variants[i].sg;
    std::cout << " Adding kmers from  " << sg.filename << std::endl;
size_t  old_size = all_kmers.size();
    KmerCompressionIndex kci(sg, max_mem*1024L*1024L*1024L);
    kci.index_graph();
    auto gk = kci.graph_kmers;
    for (auto k:gk){
       // if (k.count == 1){ - think in this sense its kmers unique to variety, not unique in graph
            all_kmers[k.kmer].insert(i);
       // }
        variants[i].kmer_total.insert(k.kmer);

    }
    size_t  new_size = all_kmers.size() - old_size;

    std::cout << "Added " << std::to_string(new_size) << " novel kmers from " << gk.size() << " kmers in assembly" << std::endl;

    return  gk.size();
};
