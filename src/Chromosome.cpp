//
// Chromosome.cpp
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


#include "Chromosome.hpp"
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <sstream>


using namespace std;


//
// HaplotypeChunk
//


ostream& operator<<(ostream& os, const HaplotypeChunk& x)
{
    os << "(" << x.position << "," << x.id << ")";
    return os;
}


istream& operator>>(istream& is, HaplotypeChunk& x)
{
    char open, comma, close;
    is >> open >> x.position >> comma >> x.id >> close;
    if (!is) return is;

    if (open != '(' ||
        comma != ',' ||
        close != ')')
        throw runtime_error("[operator>>(HaplotypeChunk)] Invalid format.");
    return is;
}


bool operator<(const HaplotypeChunk& a, const HaplotypeChunk& b)
{
    return a.position < b.position ||
           a.position == b.position && a.id < b.id;
}


bool operator>(const HaplotypeChunk& a, const HaplotypeChunk& b)
{
    return b<a;
}


bool operator==(const HaplotypeChunk& a, const HaplotypeChunk& b)
{
    return a.position == b.position && a.id == b.id;
}


bool operator!=(const HaplotypeChunk& a, const HaplotypeChunk& b)
{
    return !(a==b);
}


//
// ChromosomeEncodedID
//


ChromosomeEncodedID::ChromosomeEncodedID(unsigned int encoded)
{
    population = (encoded & 0xf0000000) >> 28;
    individual = (encoded & 0x0fffffc0) >> 6;
    pair = (encoded & 0x3e) >> 1;
    which = (encoded & 1);
}


ChromosomeEncodedID::operator unsigned int() const
{
    const unsigned int maxPopulation = 1<<4;
    const unsigned int maxIndividual = 1<<22;
    const unsigned int maxPair = 1<<5;
    const unsigned int maxWhich = 1<<1;
    
    if (population >= maxPopulation) throw runtime_error("[ChromosomeEncodedID] maxPopulation exceeded.");
    if (individual >= maxIndividual) throw runtime_error("[ChromosomeEncodedID] maxIndividual exceeded.");
    if (pair >= maxPair) throw runtime_error("[ChromosomeEncodedID] maxPair exceeded.");
    if (which >= maxWhich) throw runtime_error("[ChromosomeEncodedID] maxWhich exceeded.");

    return (population << 28) |
           (individual << 6) |
           (pair << 1) | 
           (which);
}


ostream& operator<<(ostream& os, const ChromosomeEncodedID& x)
{
    os << "<" << x.population << "," << x.individual << "," << x.pair << "," << x.which << ">";
    return os;
}


istream& operator>>(istream& is, ChromosomeEncodedID& x)
{
    char open, comma1, comma2, comma3, close;
    is >> open >> x.population >> comma1 >> x.individual >> comma2 >> x.pair >> comma3 >> x.which >> close;
    if (!is) return is;

    if (open != '<' ||
        comma1 != ',' ||
        comma2 != ',' ||
        comma3 != ',' ||
        close != '>')
        throw runtime_error("[operator>>(ChromosomeEncodedID)] Invalid format.");
    
    return is;
}


//
// Chromosome
//


Chromosome::Chromosome(unsigned int id)
{
    haplotype_chunks_.push_back(HaplotypeChunk(0, id));
}


Chromosome::Chromosome(const HaplotypeChunks& haplotype_chunks)
:   haplotype_chunks_(haplotype_chunks)
{}


Chromosome::Chromosome(const Chromosome& x, const Chromosome& y, const vector<unsigned int>& positions)
{
    bool copy_from_x = true; // false == copy from y
    size_t position_previous = 0;

    for (vector<unsigned int>::const_iterator position=positions.begin(); position!=positions.end(); ++position)
    {
        const Chromosome* p = copy_from_x ? &x : &y;
        p->extract_haplotype_chunks(position_previous, *position, haplotype_chunks_);

        copy_from_x = !copy_from_x; 
        position_previous = *position;
    }

    const Chromosome* p = copy_from_x ? &x : &y;
    p->extract_haplotype_chunks(position_previous, numeric_limits<unsigned int>::max(), haplotype_chunks_);
}


void Chromosome::extract_haplotype_chunks(unsigned int position_begin, 
                                unsigned int position_end,
                                HaplotypeChunks& result) const
{
    HaplotypeChunks::const_iterator begin = lower_bound(haplotype_chunks_.begin(), haplotype_chunks_.end(), HaplotypeChunk(position_begin,0));
    HaplotypeChunks::const_iterator end = lower_bound(haplotype_chunks_.begin(), haplotype_chunks_.end(), HaplotypeChunk(position_end,0));

    if (begin == haplotype_chunks_.end() || position_begin < begin->position)
    {
        if (begin-1 < haplotype_chunks_.begin()) throw runtime_error("[Chromosome::extract_haplotype_chunks()] Blech.");
        result.push_back(HaplotypeChunk(position_begin, (begin-1)->id));
    }

    copy(begin, end, back_inserter(result));
}


namespace
{

struct ComparePosition
{
    bool operator()(const HaplotypeChunk& a, const HaplotypeChunk& b) const {return a.position < b.position;}
};

} // namespace


HaplotypeChunks::iterator Chromosome::find_haplotype_chunk(unsigned int position, size_t index_begin)
{
    if (index_begin >= haplotype_chunks_.size()) throw runtime_error("[Chromosome::find_haplotype_chunk()] Bad index_begin.");

    HaplotypeChunks::iterator it = upper_bound(haplotype_chunks_.begin() + index_begin, haplotype_chunks_.end(), 
                                               HaplotypeChunk(position, 0), ComparePosition()); // returns first haplotype_chunk with higher position
    return it-1;
}


HaplotypeChunks::const_iterator Chromosome::find_haplotype_chunk(unsigned int position, size_t index_begin) const
{
    if (index_begin >= haplotype_chunks_.size()) throw runtime_error("[Chromosome::find_haplotype_chunk()] Bad index_begin.");

    HaplotypeChunks::const_iterator it = upper_bound(haplotype_chunks_.begin() + index_begin, haplotype_chunks_.end(), 
                                               HaplotypeChunk(position, 0), ComparePosition()); // returns first haplotype_chunk with higher position
    return it-1;
}


void Chromosome::read(istream& is)
{
    size_t haplotype_chunk_count = 0;
    is.read((char*)&haplotype_chunk_count, sizeof(size_t));
    if (haplotype_chunk_count > 10000) throw runtime_error("[Chromosome::read()] Bad haplotype_chunk_count.");
    haplotype_chunks_.resize(haplotype_chunk_count);
    is.read((char*)&haplotype_chunks_[0], sizeof(HaplotypeChunk)*haplotype_chunk_count);
}


void Chromosome::write(ostream& os) const
{
    size_t haplotype_chunk_count = haplotype_chunks_.size();
    os.write((const char*)&haplotype_chunk_count, sizeof(size_t));
    os.write((const char*)&haplotype_chunks_[0], sizeof(HaplotypeChunk)*haplotype_chunk_count);
}


bool operator==(const Chromosome& a, const Chromosome& b)
{
    if (a.haplotype_chunks().size() != b.haplotype_chunks().size()) return false;

    for (HaplotypeChunks::const_iterator it=a.haplotype_chunks().begin(), jt=b.haplotype_chunks().begin(); it!=a.haplotype_chunks().end(); ++it,++jt)
        if (*it != *jt) return false;

    return true;
}


bool operator!=(const Chromosome& a, const Chromosome& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const Chromosome& x)
{
    os << "{ ";
    copy(x.haplotype_chunks().begin(), x.haplotype_chunks().end(), ostream_iterator<HaplotypeChunk>(os, " "));
    os << "}";
    return os;
}


istream& operator>>(istream& is, Chromosome& x)
{
    string buffer;
    getline(is, buffer, '}');
    if (!is) return is;

    istringstream iss(buffer);
    char open;
    iss >> open;
    if (open != '{')
        throw runtime_error("[operator>>(Chromosome)] Invalid format.");

    HaplotypeChunks temp;
    copy(istream_iterator<HaplotypeChunk>(iss), istream_iterator<HaplotypeChunk>(), back_inserter(temp));
    x = Chromosome(temp);

    return is;
}


