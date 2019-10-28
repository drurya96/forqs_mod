//
// Population_Organisms.cpp
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


#include "Population_Organisms.hpp"
#include "Random.hpp"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iterator>
#include <set>
#include <fstream>
#include <cmath>


using namespace std;


Population_Organisms::Population_Organisms(const string& filename)
{
    if (!filename.empty())
    {
        ifstream is(filename.c_str());
        if (!is)
            throw runtime_error("[Population_Organisms] Unable to open file " + filename);
        read_text(is);
    }
}


Population_Organisms::Population_Organisms(const Organisms& organisms)
:   organisms_(organisms)
{
    if (!organisms_.empty())
    {
        chromosome_pair_count_ = organisms_.front().chromosomePairs().size();
        population_size_ = organisms_.size();

        for (Organisms::const_iterator it=organisms.begin(); it!=organisms.end(); ++it)
            if (it->chromosomePairs().size() != chromosome_pair_count_)
                throw runtime_error("[Population_Organisms] Chromosome pair counts do not match.");
    }
}


void Population_Organisms::allocate_memory()
{
    organisms_.clear();
    organisms_.resize(population_size_);

    for (Organisms::iterator it=organisms_.begin(); it!=organisms_.end(); ++it)
        it->chromosomePairs().resize(chromosome_pair_count_);
}


ChromosomePairRangeIterator Population_Organisms::begin() 
{
    return ChromosomePairRangeIterator(organisms_.begin());
}


const ChromosomePairRangeIterator Population_Organisms::begin() const
{
    // remove const on underlying container (Organisms) temporarily:
    // - this method is const, hence organisms_ is const
    // - ChromosomePairRangeIterator constructor needs non-const Organisms::iterator
    // - this method returns const ChromosomePairRangeIterator, so this is ok
    return ChromosomePairRangeIterator(const_cast<Organisms&>(organisms_).begin());
}


ChromosomePairRangeIterator Population_Organisms::end()
{
    return ChromosomePairRangeIterator(organisms_.end());
}


const ChromosomePairRangeIterator Population_Organisms::end() const
{
    // remove const on underlying container (Organisms) temporarily
    return ChromosomePairRangeIterator(const_cast<Organisms&>(organisms_).end());
}


ChromosomePairRange Population_Organisms::chromosome_pair_range(size_t organism_index)
{
    return ChromosomePairRange(organisms_[organism_index].chromosomePairs());
}


const ChromosomePairRange Population_Organisms::chromosome_pair_range(size_t organism_index) const
{
    return ChromosomePairRange(const_cast<Organism&>(organisms_[organism_index]).chromosomePairs());
}


