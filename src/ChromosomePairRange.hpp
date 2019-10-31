//
// ChromosomePairRange.hpp
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


#ifndef _CHROMOSOMEPAIRRANGE_HPP_
#define _CHROMOSOMEPAIRRANGE_HPP_


#include "Organism.hpp"


class ChromosomePairRange
{
    public:

    ChromosomePairRange(ChromosomePair* begin = 0, ChromosomePair* end = 0)
    :   begin_(begin), end_(end)
    {}

    ChromosomePairRange(ChromosomePairs& chromosome_pairs);
    
    ChromosomePair* begin() {return begin_;}
    const ChromosomePair* begin() const {return begin_;}

    ChromosomePair* end() {return end_;}
    const ChromosomePair* end() const {return end_;}

    size_t size() const {return end_ - begin_;}

    void step(int step_size) {begin_ += step_size; end_ += step_size;}


    // In-place versions of Organism constructors.

    void create_child(unsigned int id0, unsigned int id1);

    void create_child(const ChromosomePairRange& mom,
                      const ChromosomePairRange& dad,
                      const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array);

    bool equals(const ChromosomePairRange& that) const; // deep equality comparison

    private:

    ChromosomePair* begin_;
    ChromosomePair* end_;
};

class ChromosomePairRangeIterator
{
    public:

    ChromosomePairRangeIterator(Organisms::iterator it);
    ChromosomePairRangeIterator(ChromosomePair* begin, size_t chromosome_pair_count = 0);

    ChromosomePairRange& operator*();
    const ChromosomePairRange& operator*() const;

    ChromosomePairRange* operator->();
    const ChromosomePairRange* operator->() const;

    void operator++() const;

    private:

    const bool using_organism_implementation_;
    mutable Organisms::iterator current_organism_;
    mutable ChromosomePairRange current_;
    const size_t chromosome_pair_count_;

    void update_current_from_current_organism() const;
    friend bool operator==(const ChromosomePairRangeIterator& a, const ChromosomePairRangeIterator& b);
};


bool operator==(const ChromosomePairRangeIterator& a, const ChromosomePairRangeIterator& b);
bool operator!=(const ChromosomePairRangeIterator& a, const ChromosomePairRangeIterator& b);


#endif // _CHROMOSOMEPAIRRANGE_HPP_

