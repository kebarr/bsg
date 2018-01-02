//
// Created by Katie Barr (EI) on 23/11/2017.
//

#include "PhaseScaffolder.h"



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
    std::vector<std::pair<size_t, size_t> > comp_sizes;
    std::ofstream o("previously_solved_contig_names.txt");
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    sum_node_tag_mappings(sg.tags);
    for (auto component:components) {

        /*for (auto n:component) {
            if (std::find(to_check.begin(), to_check.end(), n) != to_check.end()) {
                std::cout << "Checking previously solved component \n";

                o << n << std::endl;*/

        // input for other version at:  to_look_at/subgraph_names.txt
        HaplotypeScorer hs;
        auto size = component_bps(this->sg, component);
        comp_sizes.push_back(std::make_pair(size, component.size()));
        std::cout << sg.oldnames.size() << " comp size: " << component.size() << " nodes size: " << sg.nodes.size() << std::endl;
        // oldnames
        for (auto node:sg.oldnames){

                //std::cout << " " << sg.oldnames[node-1];
           std::cout << " " << node;

        }
        for (auto node:component){

            //std::cout << " " << sg.oldnames[node-1];
            std::cout << " tag: \n";
            for (auto ts : sg.tags[node]){
                std::cout << " " << ts;
            }
            std::cout << std::endl;

        }
        std::cout << std::endl;
        // for each node, now have all tags mapping to it- ideallyl should total them
        // sum number of tags and number of mappings for each haplotype
        if (component.size() >= 6) {
            // part of the issue may be barcodes not covering the entire component
            /*    for (auto c:component){
                    std::cout << sg.oldnames[c] << " ";
                }
    std::cout << std::endl;*/
// should
            auto bubbles = sg.find_bubbles(component);
            for (auto bubble:bubbles){
                std::cout << "bubble: ";
                for (auto node:bubble) {
                    std::cout << " " << sg.oldnames[node];

                }
                std::cout << std::endl;
            }
            if (bubbles.size() > 1) {
                if (bubbles.size() <= max_bubbles) {
                    phaseable += 1;

                    int p = phase_component(bubbles, hs);
                        phased += p;

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

    std::cout << "Phased " << phased << " of " << phaseable << " phaseable components, " << too_large
              << " were too large " << " and " << not_phaseable << " did not contain enough bubbles" << std::endl;

}


void PhaseScaffolder::sum_node_tag_mappings(std::vector< std::vector<prm10xTag_t> > tag_mappings){
    sgNodeID_t counter = 0;
    for (auto n: tag_mappings){
        for (auto tag:n) {
            node_tag_mappings[counter][tag] += 1;
        }
        counter += 1;
    }
    std::cout << "summed tag mappings " << std::endl;

};

int PhaseScaffolder::phase_component(std::vector<std::vector<sgNodeID_t >> bubbles, HaplotypeScorer &hs){
    hs.find_possible_haplotypes(bubbles);
    std::cout << "mapper.reads_in_node.size()  " << mapper.reads_in_node.size() << std::endl;
    // with tags mapping to each node, just score by summing for each haplotype

    hs.count_barcode_votes(mapper);
    hs.decide_barcode_haplotype_support();
// now have mappings and barcode support
    if (hs.barcode_haplotype_mappings.size() > 0) {
        int p = hs.score_haplotypes(sg.oldnames);
        std::cout << "scored haplotypes " << p << std::endl;
        // p indicates whether scoring was successful, partially succesful or failed
        return p;
    }
}