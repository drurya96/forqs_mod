//
// Population_Organisms_Test.cpp
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
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void testPopulation_initial()
{
    if (os_) *os_ << "testPopulation_initial()\n";

    Population::Config config0;
    config0.population_size = 10;
    config0.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    Population_Organisms p0;
    p0.create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    unit_assert(p0.organisms().size() == 10);
    if (os_) *os_ << "p0:\n" << p0 << endl;

    Population::Config config1;
    config1.population_size = 4;
    config1.id_offset = 1000; // individuals==1000,1001,..., 
    config1.chromosome_pair_count = 3;
    Population_Organisms p1;
    p1.create_organisms(config1, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    unit_assert(p1.organisms().size() == 4);
    if (os_) *os_ << "p1:\n" << p1 << endl;
    if (os_) *os_ << endl;
}


void testPopulation_generated()
{
    if (os_) *os_ << "testPopulation_generated()\n";

    Population::Config config0;
    config0.population_size = 10;
    config0.chromosome_pair_count = 1;

    Population::Config config1;
    config1.population_size = 10;
    config1.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    boost::shared_ptr<Population> p0(new Population_Organisms());
    p0->create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    boost::shared_ptr<Population> p1(new Population_Organisms());
    p1->create_organisms(config1, PopulationPtrs(), PopulationDataPtrs(), rpgs);

    if (os_) *os_ << "p0:\n" << *p0 << endl;
    if (os_) *os_ << "p1:\n" << *p1 << endl;

    PopulationPtrs populations;
    populations.push_back(p0);
    populations.push_back(p1);
    
    Population::Config config_nextgen;
    config_nextgen.population_size = 100;
    config_nextgen.chromosome_pair_count = 1;
    config_nextgen.mating_distribution.push_back(MatingDistribution::Entry(.5, 0, 0));
    config_nextgen.mating_distribution.push_back(MatingDistribution::Entry(.3, 1, 0));
    config_nextgen.mating_distribution.push_back(MatingDistribution::Entry(.2, 1, 1));

    PopulationDataPtrs dummy_data;
    for (size_t i=0; i<populations.size(); ++i)
    {
        PopulationDataPtr data(new PopulationData);
        data->population_size = populations[i]->population_size();
        dummy_data.push_back(data);
    }

    Population_Organisms nextgen;
    nextgen.create_organisms(config_nextgen, populations, dummy_data, rpgs);
    if (os_) *os_ << "nextgen:\n" << nextgen << endl << flush;

    int count00 = 0;
    int count01 = 0;
    int count11 = 0;
    
    for (vector<Organism>::const_iterator it=nextgen.organisms().begin(); it!=nextgen.organisms().end(); ++it)
    {
        ChromosomeEncodedID id1(it->chromosomePairs()[0].first.haplotype_chunks()[0].id);
        ChromosomeEncodedID id2(it->chromosomePairs()[0].second.haplotype_chunks()[0].id);

        if (id1.population==0 && id2.population==0) count00++;
        else if (id1.population==1 && id2.population==1) count11++;
        else count01++;
    }

    if (os_) *os_ << "count00: " << count00 << endl;
    if (os_) *os_ << "count01: " << count01 << endl;
    if (os_) *os_ << "count11: " << count11 << endl;
    if (os_) *os_ << endl;
}


void testPopulation_fitness_constructor()
{
    if (os_) *os_ << "testPopulation_fitness_constructor()\n" << flush;

    Random::seed(0);

    // create a new population

    const size_t N = 1000;

    Population::Config config0;
    config0.population_size = 2 * N;
    config0.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    PopulationPtr p0(new Population_Organisms());
    p0->create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);

    PopulationPtrs populations;
    populations.push_back(p0);

    // fitness vector:  (1, ..., 1, 2, ..., 2)

    DataVectorPtr fitness_vector(new DataVector(config0.population_size));
    for (size_t i=0; i<N; ++i) fitness_vector->at(i) = 1;
    for (size_t i=N; i<2*N; ++i) fitness_vector->at(i) = 2;

    string id_fitness_function("id_fitness_function");
    PopulationDataPtrs population_datas;
    PopulationDataPtr data(new PopulationData);
    data->population_size = config0.population_size;
    population_datas.push_back(data);
    (*data->trait_values)[id_fitness_function] = fitness_vector;

    Population::Config config1;
    config1.population_size = config0.population_size;
    config1.chromosome_pair_count = config0.chromosome_pair_count;
    config1.mating_distribution.push_back(MatingDistribution::Entry(1, 0, 0)); // random mating
    config1.mating_distribution.default_fitness_function = id_fitness_function;

    Population_Organisms p1;
    p1.create_organisms(config1, populations, population_datas, rpgs);
    unit_assert(p1.organisms().size() == 2*N);

    size_t count1 = 0;
    size_t count2 = 0;

    for (Organisms::const_iterator it=p1.organisms().begin(); it!=p1.organisms().end(); ++it)
    {
        // cout << *it << endl; 

        unsigned int id0 = it->chromosomePairs()[0].first.haplotype_chunks()[0].id;
        unsigned int id1 = it->chromosomePairs()[0].second.haplotype_chunks()[0].id;

        if (id0 < 2*N) 
            ++count1;
        else
            ++count2;

        if (id1 < 2*N) 
            ++count1;
        else
            ++count2;
    }

    double ratio = double(count1)/count2;

    if (os_)
    {
        *os_ << "count1: " << count1 << endl
             << "count2: " << count2 << endl
             << "ratio: " << ratio << endl;
    }

    unit_assert_equal(ratio, .5, .02); // count1:count2 should be close to 1:2
}


void testPopulation_fitness_constructor_2()
{
    if (os_) *os_ << "testPopulation_fitness_constructor_2()\n";

    // create a new population

    const size_t N = 10;

    Population::Config config0;
    config0.population_size = N;
    config0.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    PopulationPtr p0(new Population_Organisms());
    p0->create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);

    if (os_) *os_ << "p0:\n" << *p0 << endl;

    PopulationPtrs populations;
    populations.push_back(p0);

    const unsigned int seed = 123;
    Random::seed(seed);

    // create new generation with fitness vector (1,...,1)

    DataVectorPtr fitness_vector(new DataVector(config0.population_size));
    for (size_t i=0; i<N; ++i) fitness_vector->at(i) = 1;

    string id_fitness_function("id_fitness_function");
    PopulationDataPtrs population_datas;
    PopulationDataPtr data(new PopulationData);
    data->population_size = config0.population_size;
    population_datas.push_back(data);
    (*data->trait_values)[id_fitness_function] = fitness_vector;

    Population::Config config1;
    config1.population_size = config0.population_size;
    config1.chromosome_pair_count = config0.chromosome_pair_count;
    config1.mating_distribution.push_back(MatingDistribution::Entry(1, 0, 0)); // random mating

    Population_Organisms p1;
    p1.create_organisms(config1, populations, population_datas, rpgs);
    unit_assert(p1.organisms().size() == N);

    if (os_) *os_ << "p1:\n" << p1 << endl;

    // create new generation, with no fitness vector specified (uniform random index)

    Random::seed(seed); // reset seed

    population_datas.front()->trait_values->clear();

    Population_Organisms p1b;
    p1b.create_organisms(config1, populations, population_datas, rpgs);
    unit_assert(p1b.organisms().size() == N);

    unit_assert(p1 == p1b); // main test: these populations should be identical

    if (os_) *os_ << "p1b:\n" << p1b << endl;
}


void test()
{
    testPopulation_initial();
    testPopulation_generated();
    testPopulation_fitness_constructor();
    testPopulation_fitness_constructor_2();
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


