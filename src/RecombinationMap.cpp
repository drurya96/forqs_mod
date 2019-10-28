//
// RecombinationMap.cpp
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


#include "RecombinationMap.hpp"
#include "Random.hpp"
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <cmath>


using namespace std;


namespace {

unsigned int factorial(unsigned int n)
{
    unsigned int result = 1; 
    for (; n>1; n--) result *= n;
    return result;
}

} // namespace


RecombinationMap::RecombinationMap(const string& filename)
{
    // read in data file
    ifstream is(filename.c_str());
    string header;
    getline(is, header);
    copy(istream_iterator<RecombinationMap::Record>(is), 
         istream_iterator<RecombinationMap::Record>(), 
         back_inserter(records_));
    if (records_.empty())
        throw runtime_error(("[RecombinationMap] Error reading file " + filename).c_str());

    // calculate Poisson distribution for number of recombination events
    // rate == cumulative geneticMap probability == expected # of events
    double rate = records_.back().geneticMap * .01; // cM * .01 = probability
    double total = 0;
    for (unsigned int i=0; i<10; i++)
    {     
        total += exp(-rate)*pow(rate, double(i))/factorial(i);
        recombinationEventDistribution_.push_back(total);
    }
}


namespace {

struct HasLowerGeneticMap
{
    bool operator()(const RecombinationMap::Record& a, const RecombinationMap::Record& b)
    {
        return a.geneticMap < b.geneticMap;
    }
};

} // namespace


unsigned int RecombinationMap::random_position()
{
    // roll randomly into the distribution, using binary search

    double max = records_.back().geneticMap;
    double roll = Random::uniform_real(0, max);

    Records::const_iterator it = lower_bound(records_.begin(), records_.end(),
                                             Record(0, 0, roll), HasLowerGeneticMap());

    if (it == records_.begin() || it == records_.end())
        throw runtime_error("[RecombinationMap::random_position()] This isn't happening.");

    // pick a position uniformly between two map positions

    unsigned int range_begin = (it-1)->position;
    unsigned int range_end = it->position - 1;
    unsigned int result = Random::uniform_integer(range_begin, range_end);

    return result;
}


vector<unsigned int> RecombinationMap::random_positions()
{
    // random number of events, according to recombinationEventDistribution_

    double roll = Random::uniform_01();
    vector<double>::const_iterator it = lower_bound(recombinationEventDistribution_.begin(),
                                                    recombinationEventDistribution_.end(),
                                                    roll);
    size_t count = it - recombinationEventDistribution_.begin();

    // pick random positions

    vector<unsigned int> result;
    for (size_t i=0; i<count; i++)
        result.push_back(random_position());
    return result;
}


istream& operator>>(istream& is, RecombinationMap::Record& r)
{
    is >> r.position >> r.combinedRate >> r.geneticMap;
    return is;
}


ostream& operator<<(ostream& os, const RecombinationMap::Record& r)
{
    os << "(" << r.position << ", " << r.combinedRate << ", " << r.geneticMap << ")";
    return os;
}


