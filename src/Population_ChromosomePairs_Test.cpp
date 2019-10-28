//
// Population_ChromosomePairs_Test.cpp
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
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_initial()
{
    if (os_) *os_ << "test_initial()\n";

    Population::Config config0;
    config0.population_size = 10;
    config0.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    Population_ChromosomePairs p0;
    p0.create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    if (os_) *os_ << "p0:\n" << p0 << endl;
    unit_assert(p0.population_size() == 10);

    Population::Config config1;
    config1.population_size = 4;
    config1.id_offset = 1000; // individuals==1000,1001,..., 
    config1.chromosome_pair_count = 3;
    Population_ChromosomePairs p1;
    p1.create_organisms(config1, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    if (os_) *os_ << "p1:\n" << p1 << endl;
    if (os_) *os_ << endl;
    unit_assert(p1.population_size() == 4);
}


void test_generated()
{
    if (os_) *os_ << "test_generated()\n";

    RecombinationPositionGenerator_Trivial recombination_position_generator("dummy_id");

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

    shared_ptr<Population> p0(new Population_ChromosomePairs());
    p0->create_organisms(config0, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    shared_ptr<Population> p1(new Population_ChromosomePairs());
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
    dummy_data.push_back(PopulationDataPtr(new PopulationData));
    dummy_data.push_back(PopulationDataPtr(new PopulationData));
    dummy_data[0]->population_size = config0.population_size;
    dummy_data[1]->population_size = config1.population_size;

    Population_ChromosomePairs nextgen;
    nextgen.create_organisms(config_nextgen, populations, dummy_data, rpgs);
    if (os_) *os_ << "nextgen:\n" << nextgen << endl << flush;
/*
    int count00 = 0;
    int count01 = 0;
    int count11 = 0;
    
    for (vector<Organism>::const_iterator it=nextgen.organisms().begin(); it!=nextgen.organisms().end(); ++it)
    {
        Chromosome::ID id1(it->chromosomePairs()[0].first.blocks()[0].id);
        Chromosome::ID id2(it->chromosomePairs()[0].second.blocks()[0].id);

        if (id1.population==0 && id2.population==0) count00++;
        else if (id1.population==1 && id2.population==1) count11++;
        else count01++;
    }

    if (os_) *os_ << "count00: " << count00 << endl;
    if (os_) *os_ << "count01: " << count01 << endl;
    if (os_) *os_ << "count11: " << count11 << endl;
    if (os_) *os_ << endl;
*/
}


void test()
{
    test_initial();
    test_generated();
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


