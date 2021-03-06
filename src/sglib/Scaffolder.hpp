//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#ifndef SG_SCAFFOLDER_HPP
#define SG_SCAFFOLDER_HPP

#include "SequenceGraph.hpp"
#include "PairedReadMapper.hpp"
#include "KmerCompressionIndex.hpp"

class Scaffolder {

public:
    ///
    /// \param _sg A SequenceGraph to scaffold
    /// \param _rms A vector of PairedReadmapper, containing the mapping of reads to _sg
    Scaffolder(SequenceGraph &_sg, std::vector<PairedReadMapper> & _rms, KmerCompressionIndex &_kci) : sg(_sg),rmappers(_rms),kci(_kci){};


    void pop_unsupported_shortbubbles();
    void find_canonical_repeats();

    ///
    /// \param source
    /// \param dest
    /// \return count of all the reads from all mappers, linking source -> dest or source -> -dest
    uint64_t count_reads_linking(sgNodeID_t source, sgNodeID_t dest);

    ///
    /// \param source
    /// \param lib the library (i.e. mapper index)
    /// \return a vecttor with dests (signed) and how many reads link to them
    std::vector<std::pair<sgNodeID_t,uint64_t>> all_read_links(sgNodeID_t source, unsigned lib);

    //2: find unsatisfied connections
    //find read support breakpoints
    //3: path finder?


    //This is a trivial strategy for the incorporation of long reads:
    //1) create a path for each long read
    //2) combine paths (OL)

    void findLongreadPaths();
    void joinLongreadPaths();
    void applyLonreadPaths();

    SequenceGraph &sg;
    std::vector<PairedReadMapper> &rmappers;
    KmerCompressionIndex &kci;



};


#endif //SG_SCAFFOLDER_HPP
