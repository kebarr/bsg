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

std::vector<std::vector<std::vector<sgNodeID_t > > > PhaseScaffolder::split_component(std::vector<std::vector<sgNodeID_t >> bubbles, int max_bubbles=12){
    // bubbles listed in order they appear in component
    std::vector<std::vector<std::vector<sgNodeID_t > >> bubbles_split;

    int to_make = bubbles.size() / (max_bubbles);
    std::cout << "splitting componenet with " << bubbles.size() << " bubbles into " << to_make << "subcomponents" << std::endl;

    if (bubbles.size()%max_bubbles == 1) {
        // can't leave part with single bubble
        int i=0;
        for (i; i < to_make-2; i ++){

            std::vector<std::vector<sgNodeID_t >> new_bubble(&bubbles[i*max_bubbles], &bubbles[(i+1)*max_bubbles]);
            bubbles_split.push_back(new_bubble);
            std::cout << new_bubble.size() << " ";
        }
        // penultimate component should have 1 less bubble so last can have 2
        i = to_make-2;
        std::vector<std::vector<sgNodeID_t >> new_bubble(&bubbles[i*max_bubbles], &bubbles[(i+1)*max_bubbles-1]);
        bubbles_split.push_back(new_bubble);
        std::cout << new_bubble.size() << " ";

        std::vector<std::vector<sgNodeID_t >> new_bubble_final(&bubbles[(i+1)*max_bubbles-1], &bubbles[bubbles.size()-1]);
        std::cout << new_bubble_final.size() << "\n";
        return bubbles_split;

    } else {
        for (int i = 0; i < to_make - 1; i++) {
            std::vector<std::vector<sgNodeID_t >> new_bubble(&bubbles[i * max_bubbles],
                                                             &bubbles[(i + 1) * max_bubbles]);
            bubbles_split.push_back(new_bubble);
        }
    }
    return bubbles_split;

}

void PhaseScaffolder::phase_components(int max_bubbles=12, int min_barcodes_mapping=10) {

//find and phase each component of gfa
    auto components = sg.connected_components();
    int too_large = 0;
    int phaseable = 0;
    int not_phaseable = 0;
    int phased = 0;
    int partial_phased = 0;
    int not_phased = 0;

    std::vector<std::pair<size_t, size_t> > comp_sizes;
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    //sum_node_tag_mappings(sg.tags);

    for (auto component:components) {

        auto size = component_bps(sg, component);
        comp_sizes.push_back(std::make_pair(size, component.size()));
        std::cout << sg.oldnames.size() << " comp size: " << component.size() << " nodes size: " << sg.nodes.size() << std::endl;

        // for each node, now have all tags mapping to it- ideallyl should total them
        // sum number of tags and number of mappings for each haplotype
        if (component.size() >= 6) {
            // part of the issue may be barcodes not covering the entire component

            auto bubbles = sg.find_bubbles(component);

            if (bubbles.size() > 1) {
                if (bubbles.size() < max_bubbles) {

                    phaseable += 1;

                    int p = phase_component(bubbles, component, min_barcodes_mapping);
                    if (p == 1) {
                        phased += 1;
                        std::cout << std::get<0>(phased_components[phased_components.size()-1].winning_pair) << std::endl;
                        std::cout << std::get<1>(phased_components[phased_components.size()-1].winning_pair) << std::endl;

                    } else if (p == 2){
                        partial_phased += 1;
                    } else if (p == 0){
                        not_phased += 1;
                    }

                } else {
                    too_large += 1;
                    // TODO output these and work out how to split up sensibly
                    auto split = split_component(bubbles);
                    for (auto sp: split){
                        int p = phase_component(sp, component, min_barcodes_mapping);
                        if (p == 1) {
                            phased += 1;

                            std::cout << phased_components[phased_components.size()-1].tags_supporting_haplotypes[std::get<0>(phased_components[phased_components.size()-1].winning_pair)].size() << std::endl;
                            std::cout << phased_components[phased_components.size()-1].tags_supporting_haplotypes[std::get<1>(phased_components[phased_components.size()-1].winning_pair)].size() << std::endl;

                        } else if (p == 2){
                            partial_phased += 1;
                        } else if (p == 0){
                            not_phased += 1;
                        }
                    }

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
        tot += c.tags_supporting_haplotypes[std::get<0>(c.winning_pair)].size();
        tot += c.tags_supporting_haplotypes[std::get<1>(c.winning_pair)].size();

    }
    double mean = tot/(phased_components.size()*2);
    for (auto c: phased_components) {
        stdev += std::pow(c.tags_supporting_haplotypes[std::get<0>(c.winning_pair)].size() - mean, 2);
        stdev += std::pow(c.tags_supporting_haplotypes[std::get<1>(c.winning_pair)].size() - mean, 2);

    }
    stdev = std::pow(stdev / (phased_components.size()*2), 0.5);
    std::cout << "Total barcodes: " << tot << " mean: " << mean << " standard deviation: " << stdev << std::endl;
}


int PhaseScaffolder::phase_component(std::vector<std::vector<sgNodeID_t >> bubbles, std::vector<sgNodeID_t > component, int min_barcodes_mapping=2){
    MappingParams mp;
    ComponentPhaser cp(sg, mapper, component, bubbles,mp);



        std::cout << "mapper.reads_in_node.size()  " << mapper.reads_in_node.size() << std::endl;
        if (cp.possible_haplotypes.size() > 1) {
            // with tags mapping to each node, just score by summing for each haplotype
            auto barcodes_map = cp.phase();
            if (barcodes_map > min_barcodes_mapping) {
// now have mappings and barcode support
                if (barcodes_map == 1) {
                    std::cout << "success:: " << cp.tags_supporting_haplotypes[std::get<0>(cp.winning_pair)].size() << std::endl;
                    std::cout << "success:: " << cp.tags_supporting_haplotypes[std::get<1>(cp.winning_pair)].size() << std::endl;

                    this->phased_components.push_back(cp);
                } else if (barcodes_map == 2) {
                    std::cout << "syccess:: " << cp.tags_supporting_haplotypes[std::get<0>(cp.winning_pair)].size() << std::endl;
                    std::cout << "syccess:: " << cp.tags_supporting_haplotypes[std::get<1>(cp.winning_pair)].size() << std::endl;

                    this->partial_phased_components.push_back(cp);
                }
                return barcodes_map;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

void PhaseScaffolder::intersect_phasings(){
    std::cout << "intersecting " << phased_components.size() << " phasings" << std::endl;
    int min_intersections = 1; // arbitrary, test
    std::vector<ComponentPhaser> phased_components_all;
    for (auto c: phased_components){
        phased_components_all.push_back(c);
    }
    for (auto c: partial_phased_components){
        phased_components_all.push_back(c);
    }
    //  for each phasing, find intersection with other phasings, help consistency by requiring x overlaps?
    std::vector<std::set<prm10xTag_t> > phasings;
    std::vector<std::set<sgNodeID_t > > phase_blocks;

    /*phasings.push_back(std::get<0>(phased_components_all[0].barcodes_supporting_winners));
    for (auto c:std::get<0>(phased_components_all[0].barcodes_supporting_winners)){
        std::cout << c << " ";
    }
    std::cout << std::endl << " first next block: ";
    phasings.push_back(std::get<1>(phased_components_all[0].barcodes_supporting_winners));
    for (auto c:std::get<1>(phased_components_all[0].barcodes_supporting_winners)){
        std::cout << c << " ";
    }*/
    std::cout << std::endl;
    //std::cout << "phasings 0 size; " << phasings[0].size() << " phasings 1 size " << phasings[1].size() << std::endl;
    int ambiguous_phasings = 0;
    std::vector<int> not_phased;
    std::vector<std::pair<std::set<prm10xTag_t>, std::set<prm10xTag_t>  >> test;
    std::set<prm10xTag_t > barcodes;
    for (auto b: barcode_node_mappings_int){
        barcodes.insert(b.first);
    }
    std::set<prm10xTag_t > b1;
    std::set<prm10xTag_t > b1_common;
    std::set<prm10xTag_t > b2;
    std::set<prm10xTag_t > b2_common;

    int cutoff = barcodes.size()/2;
    int i = 0;
    for ( auto bar:barcodes){
        if (i < 10){
            b1_common.insert(bar);
            b1.insert(bar);
        } else  if (i < cutoff){
            b1.insert(bar);
        } else if (i < cutoff + 10){
            b2_common.insert(bar);
            b2.insert(bar);
        } else {
            b2.insert(bar);
        }
        i++;
    }

            std::set<prm10xTag_t >::const_iterator iter1(b1.begin());
    std::set<prm10xTag_t >::const_iterator iter2(b2.begin());

    std::set<prm10xTag_t >::const_iterator iter11(b1_common.begin());
    std::set<prm10xTag_t >::const_iterator iter22(b2_common.begin());

    std::cout << "barcodes: " << barcodes.size() << " b1: " << b1.size() << " b2:  " << b2.size() << std::endl;
    int counter = 0;
    for (int j = 0; j < 10  ; j++) {
        std::set<prm10xTag_t> t1;
        std::set<prm10xTag_t> t2;
        for (int i = 0; i < 12; i++) {
            std::advance(iter1, 1);
            std::advance(iter2, 1);
            if (!counter % 5 == 0) {
                t1.insert(*iter1);
                t2.insert(*iter2);

            } else {

                t2.insert(*iter1);
                t1.insert(*iter2);

            }
            if (i%3==0){
                std::advance(iter11, 1);
                std::advance(iter22, 1);
                t1.insert(*iter11);
                t2.insert(*iter22);

            }
            counter += 1;


        }

        std::pair<std::set<prm10xTag_t>, std::set<prm10xTag_t>> tp = std::make_pair(t1, t2);
        test.push_back(tp);
        std::cout << " t last" << std::get<0>(test[test.size() - 1]).size() << " t1 size " << t1.size() << " t2 size " << t2.size() <<std::endl;
    }
    std::cout << " t last" << std::get<0>(test[test.size() -1]).size() <<  std::endl;

    // w
        // ant p1 abd p2 not to overkap by a lot
    //for ()
    phasings.push_back(std::get<0>(test[0]));
    phasings.push_back(std::get<1>(test[0]));

    for (int  i =1; i < test.size(); i ++){
        //auto c1 = std::get<0>(phased_components_all[i].barcodes_supporting_winners);
        //auto c2 = std::get<1>(phased_components_all[i].barcodes_supporting_winners);
        auto c1 = std::get<0>(test[i]);
        auto c2 = std::get<1>(test[i]);

        std::vector<int> best_phasing1;// = {-1, -1};// index in list of phased blocks, size of intersection
        std::vector<int> best_phasing2;// = {-1, -1};
        best_phasing1.push_back(-1);
        best_phasing1.push_back(-1);
        best_phasing2.push_back(-1);
        best_phasing2.push_back(-1);

        std::cout << "best_phasing1[0] " << best_phasing1[0] << " " << best_phasing1[1] << " " <<  best_phasing1.size()<< std::endl;
        std::cout << "best_phasing2[0] " << best_phasing2[0] << " " << best_phasing2[1] << " " <<  best_phasing2.size() << std::endl;

        std::cout << "barcides supporting hap 1 of comp " << i << " " << c1.size() << std::endl;
        for (auto c:c1){
            std::cout << c << " ";
        }
        std::cout << std::endl;
        std::cout << "barcides supporting hap 2 of comp " << i << " " << c2.size() << std::endl;

        for (auto c:c2){
            std::cout << c << " ";
        }
        std::cout << std::endl;
        for (int j =0 ; j < phasings.size(); j++){
            auto phasing = phasings[j];
            std::vector<prm10xTag_t > int1;
            std::set_intersection(phasing.begin(), phasing.end(), c1.begin(), c1.end(), std::inserter(int1, int1.begin()));
            std::vector<prm10xTag_t > int2;
            std::set_intersection(phasing.begin(), phasing.end(), c2.begin(), c2.end(), std::inserter(int2, int2.begin()));
            std::cout << "j: " << j << " i: " << i << " int1 size: " << int1.size() << " int2 size: " << int2.size()<< " best_phasing1[0]  "<< best_phasing1[0] << " best_phasing1[1]  "<< best_phasing1[1] <<  std::endl;
            std::cout << "j: " << j << " i: " << i << " int1 size: " << int1.size() << " int2 size: " << int2.size()<< " best_phasing1[0]  "<< best_phasing2[0] << " best_phasing1[1]  "<< best_phasing2[1] <<  std::endl;
            std::cout << int1.size()  << " " << best_phasing1[1] << " "  << int1.size() +  best_phasing1[1] <<  " \n";
            std::cout << int2.size()  << " " << best_phasing2[1] << " "  << int2.size() +  best_phasing2[1] << " \n";

            if ((int)int1.size() > best_phasing1[1]){
                std::cout << "new best 1: " << j << "\n";
                          best_phasing1 = {j, int1.size()}; // this loses which phasing it was but don't think that matters
            }
            if ((int)int1.size() < best_phasing1[1] ){
                std::cout << "!!!new best 1: " << int1.size()  << " " << best_phasing1[1] << "\n";
            }
            if ((int)int2.size() > best_phasing2[1]){
                std::cout << "new best 2: " << j << "\n";

                best_phasing2 = {j, int2.size()};
            }

            if ((int)int2.size() < best_phasing2[1] ){
                std::cout << "!!! new best 2: " << int2.size()  << " " << best_phasing2[1] << "\n";
            }
            int shared1 = 0;
            int shared2 = 0;

            for (auto b:phasing){
                for (auto c:c1){
                    if (c== b){
                        shared1 += 1;
                    }
                }
                for (auto c:c2){
                    if (c== b){
                        shared2 += 1;
                    }
                }

            }
            std::cout << "j: " << j << " i: " << i << " int sizes 1:  " << shared1 << ",  " << int1.size() << " 2:  "<< shared2 << ", " << int2.size() << std::endl;


        }
        std::cout << "best_phasing1[0] " << best_phasing1[0] << " " << best_phasing1[1] << std::endl;
        std::cout << "best_phasing2[0] " << best_phasing2[0] << " " << best_phasing2[1] << std::endl;

        if (best_phasing1[0] == -1 && best_phasing1[1] == -1 ){
            phasings.push_back(c1);
        } else if (best_phasing1[0] != best_phasing2[0]) {
            if (best_phasing1[1] >= min_intersections) {
                std::cout << " adding 1 to " << best_phasing1[0] << std::endl;
                for (auto c:c1) {
                    phasings[best_phasing1[0]].insert(c);
                }
            } else {
                phasings.push_back(c1);
                std::cout << " newphasing for 1 \n";

            }
            if (best_phasing2[1] >= min_intersections) {

                for (auto c:c2) {
                    phasings[best_phasing2[0]].insert(c);
                }
            } else {
                phasings.push_back(c2);
                std::cout << " newphasing for 2 \n";

            }
        } else if (best_phasing2[0] == -1 && best_phasing2[1] == -1 ){
            phasings.push_back(c2);
        } else {
            not_phased.push_back(i);
        }

    }
    std::cout << "Phased into " << phasings.size() << " blocks, " << not_phased.size() << " not phased \n" ;
}