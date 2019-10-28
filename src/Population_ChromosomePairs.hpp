//
// Population_ChromosomePairs.hpp
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


#ifndef _POPULATION_CHROMOSOMEPAIRS_HPP_
#define _POPULATION_CHROMOSOMEPAIRS_HPP_


#include "Population.hpp"


class Population_ChromosomePairs : public Population
{
    public:

    Population_ChromosomePairs(const std::string& filename = "");

    // range iteration

    virtual void allocate_memory();

    virtual ChromosomePairRangeIterator begin();
    virtual const ChromosomePairRangeIterator begin() const;

    virtual ChromosomePairRangeIterator end();
    virtual const ChromosomePairRangeIterator end() const;

    virtual ChromosomePairRange chromosome_pair_range(size_t organism_index);
    virtual const ChromosomePairRange chromosome_pair_range(size_t organism_index) const;

    private:

    ChromosomePairs chromosome_pairs_;
};


#endif // _POPULATION_CHROMOSOMEPAIRS_HPP_

