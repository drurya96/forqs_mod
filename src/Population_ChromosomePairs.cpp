//
// Population_ChromosomePairs.cpp
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


#include "Population_ChromosomePairs.hpp"
#include "Random.hpp"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iterator>
#include <set>
#include <fstream>
#include <cmath>


using namespace std;


Population_ChromosomePairs::Population_ChromosomePairs(const string& filename)
{
    if (!filename.empty())
    {
        ifstream is(filename.c_str());
        if (!is)
            throw runtime_error("[Population_ChromosomePairs] Unable to open file " + filename);
        read_text(is);
    }
}


void Population_ChromosomePairs::allocate_memory()
{
    chromosome_pairs_.clear();
    chromosome_pairs_.resize(population_size_ * chromosome_pair_count_);
}


ChromosomePairRangeIterator Population_ChromosomePairs::begin() 
{
    if (empty()) return ChromosomePairRangeIterator(0); 
    return ChromosomePairRangeIterator(&chromosome_pairs_[0], chromosome_pair_count_);
}


const ChromosomePairRangeIterator Population_ChromosomePairs::begin() const
{
    if (empty()) return ChromosomePairRangeIterator(0); 
    return ChromosomePairRangeIterator(&const_cast<ChromosomePairs&>(chromosome_pairs_)[0], chromosome_pair_count_);
}


ChromosomePairRangeIterator Population_ChromosomePairs::end()
{
    if (empty()) return ChromosomePairRangeIterator(0); 
    return ChromosomePairRangeIterator(&chromosome_pairs_[0] + chromosome_pairs_.size());
}


const ChromosomePairRangeIterator Population_ChromosomePairs::end() const
{
    if (empty()) return ChromosomePairRangeIterator(0); 
    return ChromosomePairRangeIterator(&const_cast<ChromosomePairs&>(chromosome_pairs_)[0] + chromosome_pairs_.size());
}


ChromosomePairRange Population_ChromosomePairs::chromosome_pair_range(size_t organism_index)
{
    ChromosomePair* p = &chromosome_pairs_[0] + organism_index*chromosome_pair_count_;
    return ChromosomePairRange(p, p + chromosome_pair_count_);
}


const ChromosomePairRange Population_ChromosomePairs::chromosome_pair_range(size_t organism_index) const
{
    ChromosomePair* p = &const_cast<ChromosomePairs&>(chromosome_pairs_)[0] + organism_index*chromosome_pair_count_;
    return ChromosomePairRange(p, p + chromosome_pair_count_);
}


