//
// ChromosomePairRange.cpp
//
// Created by Darren Kessner with John Novembre
//
// Copyright (c) 2013 Regents of the University of California
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// * Neither UCLA nor the names of its contributors may be used to endorse or
// promote products derived from this software without specific prior
// written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//


#include "ChromosomePairRange.hpp"
#include <stdexcept>


using namespace std;


//
// ChromosomePairRange
//


ChromosomePairRange::ChromosomePairRange(ChromosomePairs& chromosome_pairs)
:   begin_(0), end_(0)
{
    if (!chromosome_pairs.empty())
    {
        begin_ = &chromosome_pairs.front();
        end_ = begin_ + chromosome_pairs.size();
    }
}


void ChromosomePairRange::create_child(unsigned int id0, unsigned int id1)
{
    for (ChromosomePair* p=begin_; p!=end_; ++p)
    {
        p->first.haplotype_chunks().clear();
        p->first.haplotype_chunks().push_back(HaplotypeChunk(0, id0));
        p->second.haplotype_chunks().clear();
        p->second.haplotype_chunks().push_back(HaplotypeChunk(0, id1));
    }
}


void ChromosomePairRange::create_child(const ChromosomePairRange& mom,
                                       const ChromosomePairRange& dad,
                                       const RecombinationPositionGeneratorPtrs& recombination_position_generators)
{
    if (mom.size() != dad.size())
        throw runtime_error("[ChromosomePairRange::create_child()] Parents chromosome counts differ.");

    if (mom.size() != this->size())
        throw runtime_error("[ChromosomePairRange::create_child()] Parents chromosome counts differ from child.");

    if (recombination_position_generators.size() != 2)
        throw runtime_error("[ChromosomePairRange::create_child()] Recombination position generator count != 2.");

    size_t chromosome_pair_index = 0;
    const ChromosomePair* p_mom = mom.begin();
    const ChromosomePair* p_dad = dad.begin();
    ChromosomePair* p_baby = this->begin();

    for (; p_mom!=mom.end(); ++p_mom, ++p_dad, ++p_baby, ++chromosome_pair_index)
    {
        vector<unsigned int> positions_mom = recombination_position_generators[0]->get_positions(chromosome_pair_index);
        vector<unsigned int> positions_dad = recombination_position_generators[1]->get_positions(chromosome_pair_index);

        Chromosome chromosome_mom(p_mom->first, p_mom->second, positions_mom);
        p_baby->first.haplotype_chunks().swap(chromosome_mom.haplotype_chunks());

        Chromosome chromosome_dad(p_dad->first, p_dad->second, positions_dad);
        p_baby->second.haplotype_chunks().swap(chromosome_dad.haplotype_chunks());
    }
}


bool ChromosomePairRange::equals(const ChromosomePairRange& that) const
{
    if (size() != that.size()) return false;

    for (const ChromosomePair* p=begin(), *q=that.begin(); p!=end(); ++p, ++q)
        if (*p != *q) return false;

    return true;
}


//
// ChromosomePairRangeIterator
//


// Implementation note:  const is used in the interface to determine the
// result of dereferencing (ChomosomePair& vs const ChromosomePair&).
// However, the increment behavior should remain the same whether the 
// iterator is const or non-const, so the internal state is declared
// mutable.


ChromosomePairRangeIterator::ChromosomePairRangeIterator(Organisms::iterator it)
:   using_organism_implementation_(true),
    current_organism_(it), 
    chromosome_pair_count_(0)
{
    // The "organism" implementation assumes that the ChromosomePair array is held by
    // an Organism, which itself belongs to a container (Organisms). 
    // current_organism_ handles iteration through the container.
    // current_ is used as a cache, which is filled on dereference and invalidated on increment.
}


ChromosomePairRangeIterator::ChromosomePairRangeIterator(ChromosomePair* begin, 
                                                         size_t chromosome_pair_count)
:   using_organism_implementation_(false), 
    current_(begin, begin + chromosome_pair_count),
    chromosome_pair_count_(chromosome_pair_count)
{
    // The "pool" implementation assumes that the ChromosomePairs are held in a 2D array
    // where a row represents an Organism.  The chromosome_pair_count is stored as chromosome_pair_count_,
    // which is used to iterate row by row.
}


ChromosomePairRange& ChromosomePairRangeIterator::operator*() 
{
    if (using_organism_implementation_ && (current_.begin() == 0)) update_current_from_current_organism();
    return current_;
}


const ChromosomePairRange& ChromosomePairRangeIterator::operator*() const
{
    if (using_organism_implementation_ && (current_.begin() == 0)) update_current_from_current_organism();
    return current_;
}


ChromosomePairRange* ChromosomePairRangeIterator::operator->()
{
    if (using_organism_implementation_ && (current_.begin() == 0)) update_current_from_current_organism();
    return &current_;
}


const ChromosomePairRange* ChromosomePairRangeIterator::operator->() const
{
    if (using_organism_implementation_ && (current_.begin() == 0)) update_current_from_current_organism();
    return &current_;
}


void ChromosomePairRangeIterator::operator++() const
{
    if (using_organism_implementation_)
    {
        current_ = ChromosomePairRange(); // invalidate cache
        ++current_organism_;
    }
    else
    {
        if (chromosome_pair_count_ == 0)
            throw runtime_error("[ChromosomePairRangeIterator::operator++()] Step not initialized properly.");

        current_.step(chromosome_pair_count_);
    }
}


void ChromosomePairRangeIterator::update_current_from_current_organism() const
{
    if (current_organism_->chromosomePairs().empty())
        throw runtime_error("[ChromosomePairRangeIterator::update_current_from_current_organism()] No chromosome pairs.");

    // update cache
    current_ = ChromosomePairRange(current_organism_->chromosomePairs());
}


bool operator==(const ChromosomePairRangeIterator& a, const ChromosomePairRangeIterator& b)
{
    if (a.using_organism_implementation_ && b.using_organism_implementation_)
    {
        return a.current_organism_ == b.current_organism_;
    }
    else if (!a.using_organism_implementation_ && !b.using_organism_implementation_)
    {
        // Note: We only check whether beginning of range matches.
        // This facilitates range checking with automatic type conversion:
        //    ChromosomePair* -> ChromosomePairRangeIterator
        // which sets step to 0

        return a.current_.begin() == b.current_.begin();
    }

    throw runtime_error("[operator==(ChromosomePairRangeIterator)] Invalid iterator comparison: implementations differ.");
}


bool operator!=(const ChromosomePairRangeIterator& a, const ChromosomePairRangeIterator& b)
{
    return !(a==b);
}



