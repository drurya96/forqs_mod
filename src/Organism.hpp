//
// Organism.hpp
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


#ifndef _ORGANISM_HPP_
#define _ORGANISM_HPP_


#include "Chromosome.hpp"
#include "RecombinationPositionGenerator.hpp"
#include "shared_ptr.hpp"
#include <vector>


class Organism
{
    public:

    typedef std::vector<Chromosome> Gamete;

    Organism(unsigned int id0 = 0, unsigned int id1 = 0, size_t chromosomeCount = 1);

    Organism(const Gamete& g1, const Gamete& g2);

    Organism(const Organism& mom, const Organism& dad, 
             const RecombinationPositionGenerator& recombination_position_generator);

    Organism(ChromosomeEncodedID id, size_t chromosomeCount); // constructor for testing only

    const ChromosomePairs& chromosomePairs() const {return chromosomePairs_;}
    ChromosomePairs& chromosomePairs() {return chromosomePairs_;}

    Gamete create_gamete(const RecombinationPositionGenerator& recombination_position_generator) const;

    private:
    ChromosomePairs chromosomePairs_;
};


typedef std::vector<Organism> Organisms;


bool operator==(const Organism& a, const Organism& b);
bool operator!=(const Organism& a, const Organism& b);
std::ostream& operator<<(std::ostream& os, const Organism& o);
std::istream& operator>>(std::istream& is, Organism& o);


#endif // _ORGANISM_HPP_

