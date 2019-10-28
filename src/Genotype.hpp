//
// Genotype.hpp
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


#ifndef _GENOTYPER_HPP_
#define _GENOTYPER_HPP_


#include "Locus.hpp"
#include "Configurable.hpp"
#include "shared_ptr.hpp"
#include <vector>
#include <map>
#include <set>


class ChromosomePairRange;
class Organism;
class Population;
class VariantIndicator;


//
// genotype encoding:
//   (allele0, allele1) as single char
//
// e.g. genotype==0x11 means allele0==1 and allele1==1
// 


inline char genotype_make_pair(char allele_0, char allele_1)
{
    return (allele_0<<4) | allele_1;
}


inline char genotype_sum(char genotype)
{
    return (genotype>>4) + (genotype & 0x0F);
}


inline char genotype_first(char genotype)
{
    return (genotype>>4);
}


inline char genotype_second(char genotype)
{
    return (genotype & 0x0F);
}


class GenotypeData : public std::vector<char>
{
    public:

    GenotypeData() {}
    GenotypeData(const char* begin, const char* end) : std::vector<char>(begin, end) {}

    double allele_frequency() const; // note: assumes binary alleles (0/1 valued)

    // when needed: multiple allele case
    //  vector<double> multiple_allele_frequencies() const; 
    //  assume alleles are encoded as {0, ..., n-1}, return vector size n
};


typedef boost::shared_ptr<GenotypeData> GenotypeDataPtr;


class GenotypeMap : public std::map<Locus, GenotypeDataPtr> // map locus -> trait_values
{
    public:

    GenotypeDataPtr get(const Locus& locus) const // returns valid pointer, or throws
    {
        if (!count(locus))
        {
            std::ostringstream message;
            message <<  "[GenotypeMap] Locus " << locus << " not found.";
            throw std::runtime_error(message.str().c_str());
        }

        GenotypeDataPtr result = at(locus);
    
        if (!result.get())
        {
            std::ostringstream message;
            message <<  "[GenotypeMap] Null data at locus " << locus;
            throw std::runtime_error(message.str().c_str());
        }

        return at(locus);
    }
};


typedef boost::shared_ptr<GenotypeMap> GenotypeMapPtr;


//
// Genotyper
//


class Genotyper
{
    public:

    // returns genotype for a single organism at a single locus
    char genotype(const Locus& locus, 
                  const Organism& organism,
                  const VariantIndicator& indicator) const;

    // single organism, range version
    char genotype(const Locus& locus, 
                  const ChromosomePairRange& range,
                  const VariantIndicator& indicator) const;

    // genotypes population at multiple loci
    void genotype(const Loci& loci, 
                  const Population& population,
                  const VariantIndicator& indicator,
                  GenotypeMap& genotype_map) const;
};


#endif //  _GENOTYPER_HPP_

