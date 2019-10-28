//
// Chromosome.hpp
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


#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_


#include <iosfwd>
#include <vector>


struct HaplotypeChunk
{
    unsigned int position;       // position in chromosome
    unsigned int id;             // id of original DNA source

    HaplotypeChunk(unsigned int _position = 0, unsigned int _id = 0) : position(_position), id(_id) {}
};


typedef std::vector<HaplotypeChunk> HaplotypeChunks;


std::ostream& operator<<(std::ostream& os, const HaplotypeChunk& x);
std::istream& operator>>(std::istream& is, HaplotypeChunk& x);
bool operator<(const HaplotypeChunk& a, const HaplotypeChunk& b);
bool operator>(const HaplotypeChunk& a, const HaplotypeChunk& b);
bool operator==(const HaplotypeChunk& a, const HaplotypeChunk& b);
bool operator!=(const HaplotypeChunk& a, const HaplotypeChunk& b);


class Chromosome
{
    public:

    // default constructor
    Chromosome() {}

    // new chromosome, with single HaplotypeChunk
    Chromosome(unsigned int id); 

    // new chromosome with specified HaplotypeChunks, for testing
    Chromosome(const HaplotypeChunks& haplotype_chunks);

    // new chromosome via recombination
    // note:
    //  - positions must be sorted
    //  - 0 in positions <--> start with y
    Chromosome(const Chromosome& x, const Chromosome& y, const std::vector<unsigned int>& positions); 

    // access to HaplotypeChunks
    HaplotypeChunks& haplotype_chunks() {return haplotype_chunks_;}
    const HaplotypeChunks& haplotype_chunks() const {return haplotype_chunks_;}

    // append haplotype_chunks (from position_begin to position_end) to result
    void extract_haplotype_chunks(unsigned int position_begin, unsigned int position_end, HaplotypeChunks& result) const;

    // find the haplotype_chunk containing a position
    HaplotypeChunks::iterator find_haplotype_chunk(unsigned int position, size_t index_begin = 0);
    HaplotypeChunks::const_iterator find_haplotype_chunk(unsigned int position, size_t index_begin = 0) const;

    // binary read/write
    void read(std::istream& is);
    void write(std::ostream& os) const;

    private:
    HaplotypeChunks haplotype_chunks_;
};


typedef std::pair<Chromosome,Chromosome> ChromosomePair;
typedef std::vector<ChromosomePair> ChromosomePairs;




bool operator==(const Chromosome& a, const Chromosome& b);
bool operator!=(const Chromosome& a, const Chromosome& b);
std::ostream& operator<<(std::ostream& os, const Chromosome& x);
std::istream& operator>>(std::istream& is, Chromosome& x);


// Chromosome ID encoding for development/testing


struct ChromosomeEncodedID
{
    unsigned int population; // (4 bits)
    unsigned int individual; // (22 bits)
    unsigned int pair;       // (5 bits) e.g. 0-22 for humans
    unsigned int which;      // (1 bit) 
    
    ChromosomeEncodedID(unsigned int _population,
                        unsigned int _individual,
                        unsigned int _pair,
                        unsigned int _which) 
    :   population(_population),
        individual(_individual),
        pair(_pair),
        which(_which)
    {}

    // decoding from unsigned int
    explicit ChromosomeEncodedID(unsigned int encoded);

    // conversion to unsigned int
    operator unsigned int() const;
};


std::ostream& operator<<(std::ostream& os, const ChromosomeEncodedID& x);
std::istream& operator>>(std::istream& is, ChromosomeEncodedID& x);


#endif //  _CHROMOSOME_HPP_

