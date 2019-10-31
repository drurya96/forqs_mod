//
// GenotypeTest.cpp
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
#include "Population_Organisms.hpp"
#include "VariantIndicator.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


class VariantIndicator_Test : public VariantIndicator
{
    public:

    VariantIndicator_Test() : Configurable("variant_indicator_test") {}

    virtual unsigned int operator()(unsigned int chromosome_id, const Locus& locus) const
    {
        return chromosome_id;
    }
};


void test_genotype_easy()
{
    if (os_) *os_ << "test_genotype_easy()\n";

    const size_t chromosome_pair_index = 0;
    const unsigned int position = 1000000;
    Locus locus("", chromosome_pair_index, position);

    VariantIndicator_Test indicator;
    Genotyper genotyper;

    unsigned int id0 = 0;
    unsigned int id1 = 1;
    Organism mom(id0, id0, 1);
    Organism dad(id1, id1, 1);

    unsigned int genotype_mom = genotyper.genotype(locus, mom, indicator);
    unsigned int genotype_dad = genotyper.genotype(locus, dad, indicator);

    if (os_)
    {
        *os_ << "mom:\n" << mom
             << "genotype_mom: " << hex << "0x" <<  genotype_mom << dec << endl
             << "dad:\n" << dad
             << "genotype_dad: " << hex << "0x" << genotype_dad << dec << endl;
    }

    unit_assert(genotype_mom == 0);
    unit_assert(genotype_dad == 0x11);
    unit_assert(genotype_sum(genotype_dad) == 2);

    // test range version

    ChromosomePairRange range_mom(mom.chromosomePairs());
    ChromosomePairRange range_dad(dad.chromosomePairs());

    unsigned int genotype_mom_2 = genotyper.genotype(locus, range_mom, indicator);
    unsigned int genotype_dad_2 = genotyper.genotype(locus, range_dad, indicator);

    unit_assert(genotype_mom_2 == 0);
    unit_assert(genotype_dad_2 == 0x11);
    unit_assert(genotype_sum(genotype_dad) == 2);
}


void test_genotype_harder()
{
    if (os_) *os_ << "test_genotype_harder()\n";

    const size_t chromosome_pair_index = 0;
    const unsigned int position = 1000000;
    Locus locus("", chromosome_pair_index, position);

    VariantIndicator_Test indicator;
    Genotyper genotyper;

    unsigned int id0 = 0;
    unsigned int id1 = 1;

    HaplotypeChunks haplotype_chunks1;
    haplotype_chunks1.push_back(HaplotypeChunk(0, id0));
    haplotype_chunks1.push_back(HaplotypeChunk(500000, id1));
    haplotype_chunks1.push_back(HaplotypeChunk(1000000, id0)); // SNP 0 in this haplotype_chunk
    haplotype_chunks1.push_back(HaplotypeChunk(1500000, id1));
    haplotype_chunks1.push_back(HaplotypeChunk(2000000, id0));

    HaplotypeChunks haplotype_chunks2;
    haplotype_chunks2.push_back(HaplotypeChunk(0, id1));
    haplotype_chunks2.push_back(HaplotypeChunk(500000, id0));
    haplotype_chunks2.push_back(HaplotypeChunk(900000, id1)); // SNP 1 in this haplotype_chunk
    haplotype_chunks2.push_back(HaplotypeChunk(1000001, id0));
    haplotype_chunks2.push_back(HaplotypeChunk(2000000, id1));

    Chromosome chr1(haplotype_chunks1);
    Chromosome chr2(haplotype_chunks2);
    
    Organism::Gamete gamete1;
    gamete1.push_back(chr1);

    Organism::Gamete gamete2;
    gamete2.push_back(chr2);

    Organism homo1(gamete1, gamete1);
    Organism homo2(gamete2, gamete2);
    Organism hetero(gamete1, gamete2);

    unsigned int genotype_homo1 = genotyper.genotype(locus, homo1, indicator);
    unsigned int genotype_homo2 = genotyper.genotype(locus, homo2, indicator);
    unsigned int genotype_hetero = genotyper.genotype(locus, hetero, indicator);

    if (os_)
    {
        *os_ << "homo1:\n" << homo1
             << "genotype_homo1: " << hex << "0x" << genotype_homo1 << dec << endl
             << "homo2:\n" << homo2
             << "genotype_homo2: " << hex << "0x" << genotype_homo2 << dec << endl
             << "hetero:\n" << hetero
             << "genotype_hetero: " << hex << "0x" << genotype_hetero << dec << endl;
    }

    unit_assert(genotype_homo1 == 0);
    unit_assert(genotype_homo2 == 0x11);
    unit_assert(genotype_sum(genotype_homo2) == 2);
    unit_assert(genotype_hetero == 0x01);
    unit_assert(genotype_first(genotype_hetero) == 0);
    unit_assert(genotype_second(genotype_hetero) == 1);
    unit_assert(genotype_sum(genotype_hetero) == 1);

    // test genotype(population)

    Organisms organisms;
    organisms.push_back(homo1);
    organisms.push_back(homo2);
    organisms.push_back(hetero);
    organisms.push_back(homo1);
    organisms.push_back(homo2);
    organisms.push_back(hetero);

    Loci loci;
    loci.insert(locus);
    
    Population_Organisms population(organisms);

    GenotypeMap genotype_map;
    genotyper.genotype(loci, population, indicator, genotype_map);

    GenotypeDataPtr genotypes = genotype_map.at(locus);
    
    if (os_)
    {
        *os_ << "genotypes: " << genotypes->size() << endl;
        *os_ << hex;
        copy(genotypes->begin(), genotypes->end(), ostream_iterator<double>(*os_, " "));
        *os_ << dec << endl;
    }

    unit_assert(genotypes->size() == 6);
    unit_assert(genotypes->at(0) == 0 && genotypes->at(3) == 0);
    unit_assert(genotypes->at(1) == 0x11 && genotypes->at(4) == 0x11);
    unit_assert(genotypes->at(2) == 0x1 && genotypes->at(5) == 0x1);
}


void test_allele_frequency()
{
    GenotypeData data;

    // Hardy-Weinberg with p=.2, N=100
    for (size_t i=0; i<64; i++) data.push_back(0);
    for (size_t i=0; i<32; i++) data.push_back(1);
    for (size_t i=0; i<4; i++) data.push_back(2);

    const double epsilon = 1e-12;
    unit_assert_equal(data.allele_frequency(), .2, epsilon);
}


void test_map_get()
{
    GenotypeDataPtr data(new GenotypeData);
    Locus locus("locus", 1, 100000);
    Locus locus_bad("locus_bad", 2, 200000);
    GenotypeMap gm;
    gm[locus] = data;

    GenotypeDataPtr retrieved = gm.get(locus); // ok
    unit_assert(retrieved.get() == data.get());

    // try to retrieve data that isn't in the map

    bool caught = false;
    try
    {
        GenotypeDataPtr dummy = gm.get(locus_bad);
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught exception:\n" << e.what() << endl;
        caught = true;
    }
    unit_assert(caught);

    // try to retrieve null data vector

    GenotypeDataPtr null;
    gm[locus_bad] = null;
    caught = false;
    try
    {
        GenotypeDataPtr dummy = gm.get(locus_bad);
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught exception:\n" << e.what() << endl;
        caught = true;
    }
    unit_assert(caught);
}


void test()
{
    test_genotype_easy();
    test_genotype_harder();
    test_allele_frequency();
    test_map_get();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        test();
        return 0;
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


