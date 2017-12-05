//
// Created by Katie Barr (EI) on 23/11/2017.
//

#include "PhaseScaffolder.h"



PhaseScaffolder::PhaseScaffolder(SequenceGraph & sg): sg(sg), mapper(sg){
}

void PhaseScaffolder::load_mappings(std::string r1_filename, std::string r2_filename, std::string fasta_filename, uint64_t max_mem_gb, std::string to_map=""){

    mapper.map_reads(r1_filename, r2_filename, fasta_filename, prm10x, max_mem_gb);
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
    if (to_map.size()>0){
        std::cout << "Writting mappings to disk" << std::endl;
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

void PhaseScaffolder::phase_components() {

//find and phase each component of gfa
    auto components = sg.connected_components();
    int too_large = 0;
    int phaseable = 0;
    int not_phaseable = 0;
    int phased = 0;
    std::vector<sgNodeID_t > to_check = {sg.oldnames_to_ids["448660"], sg.oldnames_to_ids["461782"], sg.oldnames_to_ids["468510"], sg.oldnames_to_ids["510080"],
                                         sg.oldnames_to_ids["514764"], sg.oldnames_to_ids["843062"], sg.oldnames_to_ids["437632"], sg.oldnames_to_ids["440422"],
                                         sg.oldnames_to_ids["448660"], sg.oldnames_to_ids["437440"], sg.oldnames_to_ids["417378"], sg.oldnames_to_ids["413884"],
                                         sg.oldnames_to_ids["411428"], sg.oldnames_to_ids["402916"]};
    std::vector<std::pair <size_t , size_t > > comp_sizes;
    std::ofstream o("previously_solved_contig_names.txt");
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    for (auto component:components) {

            /*for (auto n:component) {
                if (std::find(to_check.begin(), to_check.end(), n) != to_check.end()) {
                    std::cout << "Checking previously solved component \n";

                    o << n << std::endl;*/

            // input for other version at:  to_look_at/subgraph_names.txt
            HaplotypeScorer hs;
            auto size = component_bps(this->sg, component);
            comp_sizes.push_back(std::make_pair(size, component.size()));
            if (component.size() >= 6) {
                // part of the issue may be barcodes not covering the entire component
                /*    for (auto c:component){
                        std::cout << sg.oldnames[c] << " ";
                    }
        std::cout << std::endl;*/
// should
                auto bubbles = sg.find_bubbles(component);
                if (bubbles.size() > 1) {
                    if (bubbles.size() <= 20) {
                        phaseable += 1;
                        hs.find_possible_haplotypes(bubbles);
                        std::cout << "mapper.reads_in_node.size()  " << mapper.reads_in_node.size() << std::endl;

                        hs.count_barcode_votes(mapper);
                        hs.decide_barcode_haplotype_support();
// now have mappings and barcode support
                        if (hs.barcode_haplotype_mappings.size() > 0) {
                            int p = hs.score_haplotypes(sg.oldnames);
                            std::cout << "scored haplotypes " << p << std::endl;
                            phased += p;
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

    std::cout << "Phased " << phased << " of " << phaseable << " phaseable components, " << too_large << " were too large " << " and " << not_phaseable << " did not contain enough bubbles" << std::endl;
    std::cout << "component sizes \n";
              int i = 0;
    for (auto c:comp_sizes){
        std::cout << std::get<0>(c) << ": " << std::get<1>(c) << "  ";
        i++;
        if (i%100 == 0){
            std::cout << std::endl;
        }
    }
}