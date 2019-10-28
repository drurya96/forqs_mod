//
// Organism.cpp
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


#include "Organism.hpp"
#include "Random.hpp"
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iterator>


using namespace std;


Organism::Organism(unsigned int id0, unsigned int id1, size_t chromosomeCount)
{
	cout << "AD Controlled Crash" << endl;
	cout << "Attempted to make organism from ids" << endl;
	exit(1);
    for (size_t i=0; i<chromosomeCount; ++i)
        chromosomePairs_.push_back(make_pair(Chromosome(id0), Chromosome(id1)));
}


Organism::Organism(const Gamete& g1, const Gamete& g2)
{
	cout << "AD Controlled Crash" << endl;
	cout << "Attempted to make organism from gametes" << endl;
	exit(1);
    if (g1.size() != g2.size())
        throw runtime_error("[Organism::Organism()] Different gamete sizes.");

    for (size_t i=0; i<g1.size(); ++i)
        chromosomePairs_.push_back(make_pair(g1[i],g2[i]));
}


Organism::Organism(const Organism& mom, const Organism& dad,
                   const RecombinationPositionGenerator& recombination_position_generator)
{
	cout << "AD Controlled Crash" << endl;
	cout << "Attempted to make organism from parents" << endl;
	exit(1);
    if (mom.chromosomePairs_.size() != dad.chromosomePairs_.size())
        throw runtime_error("[Organism::Organism(mom, dad)] Parents chromosome counts differ.");

    this->chromosomePairs_.reserve(mom.chromosomePairs_.size());

    size_t chromosome_index = 0;
    for (ChromosomePairs::const_iterator it=mom.chromosomePairs_.begin(), jt=dad.chromosomePairs_.begin();
         it!=mom.chromosomePairs_.end(); ++it, ++jt, ++chromosome_index)
    {
        vector<unsigned int> positions_mom = recombination_position_generator.get_positions(chromosome_index);
        vector<unsigned int> positions_dad = recombination_position_generator.get_positions(chromosome_index);

        this->chromosomePairs_.push_back(make_pair(
            Chromosome(it->first, it->second, positions_mom),
            Chromosome(jt->first, jt->second, positions_dad)));
    }
}


Organism::Organism(ChromosomeEncodedID id, size_t chromosomeCount)
{
	cout << "AD Controlled Crash" << endl;
	cout << "Attempted to make organism from random" << endl;
	exit(1);
    cerr << "[Organism] Warning: constructor for testing only\n";

    ChromosomeEncodedID id0(id);
    ChromosomeEncodedID id1(id);
    id0.which = 0;
    id1.which = 1;

    for (size_t i=0; i<chromosomeCount; ++i)
    {
        id0.pair = id1.pair = i;
        chromosomePairs_.push_back(make_pair(Chromosome(id0), Chromosome(id1)));
    }
}


Organism::Gamete Organism::create_gamete(const RecombinationPositionGenerator& recombination_position_generator) const
{

	cout << "AD Controlled Crash" << endl;
	cout << "Attempted to make gamete as organism" << endl;
	exit(1);

    Gamete result;

    for (ChromosomePairs::const_iterator it=chromosomePairs_.begin(); it!=chromosomePairs_.end(); ++it)
    {
        vector<unsigned int> positions = 
            recombination_position_generator.get_positions(it-chromosomePairs_.begin());

        result.push_back(Chromosome(it->first, it->second, positions));
    }

    return result;
}


bool operator==(const Organism& a, const Organism& b)
{
    if (a.chromosomePairs().size() != b.chromosomePairs().size())
        return false;

    for (ChromosomePairs::const_iterator it=a.chromosomePairs().begin(), jt=b.chromosomePairs().begin();
         it!=a.chromosomePairs().end(); ++it, ++jt)
        if (*it!=*jt) return false;

    return true;
}


bool operator!=(const Organism& a, const Organism& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Organism& o)
{
    if (o.chromosomePairs().empty()) throw runtime_error("[Organism::operator<<] No chromosome pairs.");

    for (ChromosomePairs::const_iterator it=o.chromosomePairs().begin(); it!=o.chromosomePairs().end(); ++it)
        os << "+ " << it->first << "\n- " << it->second << endl;
    return os;
}


istream& operator>>(istream& is, Organism& o)
{
    Organism::Gamete g_plus, g_minus;

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!is || buffer.empty()) break;

        char plus_minus;
        Chromosome chromosome(0);
        istringstream iss(buffer);
        iss >> plus_minus >> chromosome;

        if (plus_minus == '+')
            g_plus.push_back(chromosome); 
        else if (plus_minus == '-')
            g_minus.push_back(chromosome); 
        else
            throw runtime_error("[operator>>(Organism)] Invalid format.");
    }

    if (!g_plus.empty())
        o = Organism(g_plus, g_minus);

    return is;
}


