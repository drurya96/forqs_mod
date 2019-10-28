//
// Genotype.cpp
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


#include "Genotype.hpp"
#include "VariantIndicator.hpp"
#include "boost/lexical_cast.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>


using namespace std; 


//
// GenotypeData
//


double GenotypeData::allele_frequency() const
{
    vector<char> sums(size());
    transform(begin(), end(), sums.begin(), genotype_sum);
    double sum = accumulate(sums.begin(), sums.end(), 0.0);
    return sum/size()/2;
}


//
// Genotyper
//


char Genotyper::genotype(const Locus& locus, 
                         const Organism& organism,
                         const VariantIndicator& indicator) const
{
    if (locus.chromosome_pair_index >= organism.chromosomePairs().size())
        throw runtime_error("[Genotyper::genotype(Organism) chromosome_pair_index out of range.");
    const ChromosomePair& cp = organism.chromosomePairs()[locus.chromosome_pair_index];
    const HaplotypeChunk& haplotype_chunk0 = *cp.first.find_haplotype_chunk(locus.position);
    const HaplotypeChunk& haplotype_chunk1 = *cp.second.find_haplotype_chunk(locus.position);
    return genotype_make_pair(indicator(haplotype_chunk0.id, locus), indicator(haplotype_chunk1.id, locus));
}


char Genotyper::genotype(const Locus& locus, 
                         const ChromosomePairRange& range,
                         const VariantIndicator& indicator) const
{
    if (locus.chromosome_pair_index >= range.size())
        throw runtime_error("[Genotyper::genotype(ChromosomePairRange) chromosome_pair_index out of range.");
    const ChromosomePair& cp = range.begin()[locus.chromosome_pair_index];
    const HaplotypeChunk& haplotype_chunk0 = *cp.first.find_haplotype_chunk(locus.position);
    const HaplotypeChunk& haplotype_chunk1 = *cp.second.find_haplotype_chunk(locus.position);
    return genotype_make_pair(indicator(haplotype_chunk0.id, locus), indicator(haplotype_chunk1.id, locus));
}


void Genotyper::genotype(const Loci& loci, 
                         const Population& population,
                         const VariantIndicator& indicator,
                         GenotypeMap& genotype_map) const
{
    for (Loci::const_iterator locus=loci.begin(); locus!=loci.end(); ++locus)
    {
        GenotypeDataPtr genotypes(new GenotypeData);
        genotypes->reserve(population.population_size());

        const ChromosomePairRangeIterator range_end = population.end();
        for (const ChromosomePairRangeIterator range=population.begin(); range!=range_end; ++range)
            genotypes->push_back(genotype(*locus, *range, indicator));

        genotype_map[*locus] = genotypes;
    }
}


