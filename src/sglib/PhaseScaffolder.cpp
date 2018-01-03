//
// Created by Katie Barr (EI) on 23/11/2017.
//

#include "PhaseScaffolder.h"

/*
class PhasedComponent{
public:
    std::vector<sgNodeID_t> component;
    std::vector<>
};*/


PhaseScaffolder::PhaseScaffolder(SequenceGraph & sg): sg(sg), mapper(sg){
}

void PhaseScaffolder::load_mappings(std::string r1_filename, std::string r2_filename, std::string fasta_filename, uint64_t max_mem_gb, std::string to_map=""){

    mapper.map_reads(r1_filename, r2_filename, prm10x, max_mem_gb);
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
    if (to_map.size()>0){
        std::cout << "Writing mappings to disk" << std::endl;
        mapper.save_to_disk(to_map);
    }
    mapper.print_stats();
}

void PhaseScaffolder::load_mappings_from_file(std::string filename){
    std::cout << "Loading mappings from " << filename << std::endl;
    mapper.load_from_disk(filename);
    mapper.print_stats();
}


void PhaseScaffolder::output_bubbles(std::string bubble_filename) {

        std::ofstream out2(bubble_filename);
//find each component of gfa
    std::cout << "Finding components" << std::endl;
    auto components = sg.connected_components();
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    int counter = 0;
    int counter2 = 0;
    for (auto component:components) {
        if (component.size() >= 6) {
            auto bubbles = sg.find_bubbles(component);

            if (bubbles.size() > 1) {
                counter2+=1;
                for (auto bubble:bubbles) {

                    for (auto bubble_c:bubble) {
                        out2 << ">" << sg.oldnames[bubble_c] << "_" << counter << std::endl
                             << sg.nodes[bubble_c].sequence << std::endl;
                    }

                }

            }

            counter += 1;
        }
    }
    std::cout << counter << " components with enough contigs to be phaseable, " << counter2 << " contained more tan 1 bubble, i.e. are phaseable \n";
}

size_t  component_bps(SequenceGraph &sg, std::vector<sgNodeID_t > component){
    size_t res = 0;
    for (auto n:component){
        res += sg.nodes[n].sequence.size();
    }
    return  res;
}

void PhaseScaffolder::phase_components(int max_bubbles=12) {

//find and phase each component of gfa
    auto components = sg.connected_components();
    int too_large = 0;
    int phaseable = 0;
    int not_phaseable = 0;
    int phased = 0;
    int partial_phased = 0;
    int not_phased = 0;

    std::vector<std::pair<size_t, size_t> > comp_sizes;
    std::ofstream o("previously_solved_contig_names.txt");
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    sum_node_tag_mappings(sg.tags);

    for (auto component:components) {

        auto size = component_bps(this->sg, component);
        comp_sizes.push_back(std::make_pair(size, component.size()));
        std::cout << sg.oldnames.size() << " comp size: " << component.size() << " nodes size: " << sg.nodes.size() << std::endl;

        // for each node, now have all tags mapping to it- ideallyl should total them
        // sum number of tags and number of mappings for each haplotype
        if (component.size() >= 6) {
            // part of the issue may be barcodes not covering the entire component

            auto bubbles = sg.find_bubbles(component);

            if (bubbles.size() > 1) {
                if (bubbles.size() <= max_bubbles) {

                    phaseable += 1;

                    int p = phase_component(bubbles);
                    if (p == 1) {
                        phased += 1;
                    } else if (p == 2){
                        partial_phased += 1;
                    } else if (p == 0){
                        not_phased += 1;
                    }

                } else {
                    too_large += 1;
                    // TODO output these and work out how to split up sensibly

                }
            } else {
                not_phaseable += 1;
            }
        }
    }
    //  }
    //}

    std::cout << "Phased " << phased <<  " partial phased " << partial_phased << " not Phased " << not_phased <<  " of " << phaseable << " phaseable components, " << too_large
              << " were too large " << " and " << not_phaseable << " did not contain enough bubbles" << std::endl;
    intersect_phasings();
    std::cout << "Phased components: " << phased_components.size() << " partyially phased: " << partial_phased_components.size() << std::endl;
    print_barcode_stats();
}

void PhaseScaffolder::print_barcode_stats(){
    int tot = 0;
    double stdev = 0;

    for (auto c:phased_components){
        tot += std::get<0>(c.barcodes_supporting_winners).size();
        tot += std::get<1>(c.barcodes_supporting_winners).size();

    }
    double mean = tot/(phased_components.size()*2);
    for (auto c: phased_components) {
        stdev += std::pow(std::get<0>(c.barcodes_supporting_winners).size() - mean, 2);
        stdev += std::pow(std::get<1>(c.barcodes_supporting_winners).size() - mean, 2);

    }
    stdev = std::pow(stdev / (phased_components.size()*2), 0.5);
    std::cout << "Total barcodes: " << tot << " mean: " << mean << " standard deviation: " << stdev << std::endl;
}


void PhaseScaffolder::sum_node_tag_mappings(std::vector< std::vector<prm10xTag_t> > tag_mappings, int min_tag_count=1){
    std::map<sgNodeID_t, std::map<prm10xTag_t, int > > node_tag_mappings_int;
    sgNodeID_t counter = 0;
    for (auto n: tag_mappings){
        for (auto tag:n) {
            node_tag_mappings_int[counter][tag] += 1;
        }
        counter += 1;
    }
    std::cout << "summed tag mappings: " << node_tag_mappings_int.size() << std::endl;
    int discarded_barodes = 0;
    int kept_barodes = 0;
    for (auto b:node_tag_mappings_int){
        for (auto t:b.second){
            if (t.second > min_tag_count){
                node_tag_mappings[b.first][t.first] = t.second;
                kept_barodes += 1;
            } else {
                discarded_barodes += 1;
            }

        }
    }
    std::cout << "filtered tag mappings: " << node_tag_mappings.size() << " discarded " << discarded_barodes << " barcodes and kept " << kept_barodes << std::endl;

};

int PhaseScaffolder::phase_component(std::vector<std::vector<sgNodeID_t >> bubbles){
    HaplotypeScorer hs;
    std::map<sgNodeID_t, std::map<prm10xTag_t, int > > relevant_mappings;
    std::map<prm10xTag_t, std::vector<sgNodeID_t > > barcode_node_mappings;
    for (auto bubble:bubbles){
        for (auto b: bubble) {
            relevant_mappings[b] = node_tag_mappings[b];
            for (auto t: node_tag_mappings[b]){
                barcode_node_mappings[t.first].push_back(b);
            }

        }
    }
    hs.find_possible_haplotypes(bubbles);
    std::cout << "mapper.reads_in_node.size()  " << mapper.reads_in_node.size() << std::endl;
    // with tags mapping to each node, just score by summing for each haplotype
    hs.decide_barcode_haplotype_support(relevant_mappings, barcode_node_mappings);
    int res = hs.score_haplotypes(sg.oldnames);
// now have mappings and barcode support
    if (hs.barcode_haplotype_mappings.size() > 0) {
        //std::cout << "scored haplotypes " << p << std::endl;
        // p indicates whether scoring was successful, partially succesful or failed
        //return p;
    }

    if (res == 1){
        phased_components.push_back(hs);
    } else if (res == 2) {
        partial_phased_components.push_back(hs);
    }
    return res;
}

void PhaseScaffolder::intersect_phasings(){
    std::cout << "intersecting " << phased_components.size() << " phasings" << std::endl;
    int min_intersections = 10; // arbitrary, test
    //  for each phasing, find intersection with other phasings, help consistency by requiring x overlaps?
    std::vector<std::set<prm10xTag_t> > phasings;
    phasings.push_back(std::get<0>(phased_components[0].barcodes_supporting_winners));
    phasings.push_back(std::get<1>(phased_components[0].barcodes_supporting_winners));
    int ambiguous_phasings = 0;
    std::vector<int> not_phased;
    for (int  i =1; i < phased_components.size(); i ++){
        auto c1 = std::get<0>(phased_components[i].barcodes_supporting_winners);
        auto c2 = std::get<1>(phased_components[i].barcodes_supporting_winners);
        std::vector<unsigned long> best_phasing1;
        std::vector<unsigned long > best_phasing2;

        for (int j =0 ; j < phasings.size(); j++){
            auto phasing = phasings[j];
            std::vector<prm10xTag_t > int1;
            std::set_intersection(phasing.begin(), phasing.end(), c1.begin(), c1.end(), std::back_inserter(int1));
            std::vector<prm10xTag_t > int2;
            std::set_intersection(phasing.begin(), phasing.end(), c2.begin(), c2.end(), std::back_inserter(int2));
            if (int1.size() > best_phasing1[1]){
                best_phasing1 = {j, int1.size()}; // this loses which phasing it was but don't think that matters
            }
            if (int2.size() > best_phasing2[1]){
                best_phasing2 = {j, int2.size()};
            }

        }
        if (best_phasing1[0] != best_phasing2[0]) {
            if (best_phasing1[1] >= min_intersections) {

                for (auto c:c1) {
                    phasings[best_phasing1[0]].insert(c);
                }
            } else {
                phasings.push_back(c1);

            }
            if (best_phasing2[1] >= min_intersections) {

                for (auto c:c2) {
                    phasings[best_phasing2[0]].insert(c);
                }
            } else {
                phasings.push_back(c2);
            }
        } else {
            not_phased.push_back(i);
        }

    }
    std::cout << "Phased into " << phasings.size() << " blocks, " << not_phased.size() << " not phased \n" ;
}