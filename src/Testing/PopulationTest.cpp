//
// PopulationTest.cpp
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


#include "Population.hpp"
#include "Population_Organisms.hpp"
#include "Population_ChromosomePairs.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_MatingDistribution()
{
    if (os_) *os_ << "testMatingDistribution()\n" << flush;

    MatingDistribution md;
    md.push_back(MatingDistribution::Entry(.6, 0, 10));
    md.push_back(MatingDistribution::Entry(.2, 1, "ff1", 11, "ff2"));
    md.push_back(MatingDistribution::Entry(.1, 2, "ff2", 12, "ff1"));

    unit_assert(md.entries().size() == 3);
    unit_assert(md.cumulative_weights()[0] == .6);
    unit_assert(md.cumulative_weights()[1] == .8);
    unit_assert(md.cumulative_weights()[2] == .9);

    unit_assert(md.entries()[0].first == 0);
    unit_assert(md.entries()[0].first_fitness == "");
    unit_assert(md.entries()[0].second == 10);
    unit_assert(md.entries()[0].second_fitness == "");

    unit_assert(md.entries()[1].first == 1);
    unit_assert(md.entries()[1].first_fitness == "ff1");
    unit_assert(md.entries()[1].second == 11);
    unit_assert(md.entries()[1].second_fitness == "ff2");

    unit_assert(md.entries()[2].first == 2);
    unit_assert(md.entries()[2].first_fitness == "ff2");
    unit_assert(md.entries()[2].second == 12);
    unit_assert(md.entries()[2].second_fitness == "ff1");

    vector<int> counts(3);

    for (int i=0; i<9000; i++)
    {
        const MatingDistribution::Entry& entry = md.random();
        counts[entry.first]++;
    }

    for (size_t i=0; i<3; i++)
        if (os_) *os_ << "count " << i << ": " << counts[i] << endl;

    // test I/O

    ostringstream oss;
    oss << md;
    if (os_) *os_ << oss.str() << endl;

    MatingDistribution md2;
    md2.push_back(MatingDistribution::Entry(.42, 4, 5));
    if (os_) *os_ << "md2 before: " << md2 << endl;

    istringstream iss(oss.str());
    iss >> md2;
    if (os_) *os_ << "md2 after: " << md2 << endl;
    unit_assert(md==md2);

    if (os_) *os_ << endl;
}


void test_MatingDistribution_validate_entries()
{
    if (os_) *os_ << "test_MatingDistribution_validate_entries()\n" << flush;

    MatingDistribution md;
    md.default_fitness_function = "default";
    md.push_back(MatingDistribution::Entry(1, 0, 0));
    md.push_back(MatingDistribution::Entry(1, 1, "ff1", 1, "ff2"));

    DataVectorPtr zero(new DataVector(1,0));
    DataVectorPtr nonzero(new DataVector(1,1));
    unit_assert(zero->all_zero());
    unit_assert(!nonzero->all_zero());

    PopulationDataPtrs population_datas;
    population_datas.push_back(PopulationDataPtr(new PopulationData));
    population_datas.push_back(PopulationDataPtr(new PopulationData));

    // everything valid

    (*population_datas[0]->trait_values)["default"] = nonzero;
    (*population_datas[1]->trait_values)["ff1"] = nonzero;
    (*population_datas[1]->trait_values)["ff2"] = nonzero;

    md.validate_entries(population_datas);
    vector<size_t> counts(2);
    for (size_t i=0; i<10; ++i)
        ++counts[md.random().first];
    unit_assert(counts[0] && counts[1]);

    // nothing valid

    if (os_) *os_ << "everything invalid:\n";
    (*population_datas[0]->trait_values)["default"] = zero;
    (*population_datas[1]->trait_values)["ff1"] = zero;
    (*population_datas[1]->trait_values)["ff2"] = zero;
    bool caught = false;
    try
    {
        md.validate_entries(population_datas);
    }
    catch (...)
    {
        caught = true;
    }
    unit_assert(caught);

    // default not valid

    if (os_) *os_ << "default invalid:\n";
    (*population_datas[0]->trait_values)["default"] = zero;
    (*population_datas[1]->trait_values)["ff1"] = nonzero;
    (*population_datas[1]->trait_values)["ff2"] = nonzero;
    md.validate_entries(population_datas);
    counts = vector<size_t>(2); // reset
    for (size_t i=0; i<10; ++i)
        ++counts[md.random().first];
    unit_assert(!counts[0] && counts[1]);

    // ff1 not valid

    if (os_) *os_ << "ff1 invalid:\n";
    (*population_datas[0]->trait_values)["default"] = nonzero;
    (*population_datas[1]->trait_values)["ff1"] = zero;
    (*population_datas[1]->trait_values)["ff2"] = nonzero;
    md.validate_entries(population_datas);
    counts = vector<size_t>(2); // reset
    for (size_t i=0; i<10; ++i)
        ++counts[md.random().first];
    unit_assert(counts[0] && !counts[1]);

    // ff2 not valid

    if (os_) *os_ << "ff2 invalid:\n";
    (*population_datas[0]->trait_values)["default"] = nonzero;
    (*population_datas[0]->trait_values)["ff1"] = nonzero;
    (*population_datas[0]->trait_values)["ff2"] = zero;
    md.validate_entries(population_datas);
    counts = vector<size_t>(2); // reset
    for (size_t i=0; i<10; ++i)
        ++counts[md.random().first];
    unit_assert(counts[0] && !counts[1]);
}


void test_Population_Config_IO()
{
    if (os_) *os_ << "test_Population_Config_IO()\n";

    Population::Config config1;
    config1.population_size = 4;
    config1.id_offset = 1000; // individuals==1000,1001,..., 
    config1.chromosome_pair_count = 3;
    config1.mating_distribution.push_back(MatingDistribution::Entry(.42, 0,0));
    config1.mating_distribution.push_back(MatingDistribution::Entry(.66, 1,0));
    config1.mating_distribution.push_back(MatingDistribution::Entry(.23, 2,1));

    ostringstream oss;
    oss << config1;

    if (os_) *os_ << "config1:\n" << config1 << endl;

    Population::Config config2;
    unit_assert(config1 != config2);
    istringstream iss(oss.str());
    iss >> config2;

    if (os_) *os_ << "config2:\n" << config2 << endl;
    unit_assert(config1 == config2);

    if (os_) *os_ << endl;
}


void test_generations_IO()
{
    if (os_) *os_ << "test_generation_configs_IO()\n" << flush;

    vector<Population::Configs> populationConfigs;

    const size_t chromosome_pair_count_ = 3;
    const unsigned int populationSize_ = 10000;
    const double admixtureProportion_ = .8; // fraction of genes from 1st population

    // generation 0 (ancestral populations)

    populationConfigs.push_back(vector<Population::Config>(3));

    Population::Config* config_pop = &populationConfigs[0][0];
    config_pop->population_size = 0;

    config_pop = &populationConfigs[0][1];
    config_pop->population_size = populationSize_;
    config_pop->chromosome_pair_count = chromosome_pair_count_;

    config_pop = &populationConfigs[0][2];
    config_pop->population_size = populationSize_;
    config_pop->chromosome_pair_count = chromosome_pair_count_;

    // generation 1 (initial admixture)

    populationConfigs.push_back(vector<Population::Config>(1));

    config_pop = &populationConfigs[1][0];
    config_pop->population_size = populationSize_;
    double p = admixtureProportion_;
    config_pop->mating_distribution.push_back(MatingDistribution::Entry(p*p, 1,1));
    config_pop->mating_distribution.push_back(MatingDistribution::Entry(2*p*(1-p), 1,2));
    config_pop->mating_distribution.push_back(MatingDistribution::Entry((1-p)*(1-p), 2,2));

    // subsequent generations - just recombination

    const size_t generation_count = 8;

    for (size_t generation=2; generation<generation_count; generation++)
    {
        populationConfigs.push_back(vector<Population::Config>(1));
        config_pop = &populationConfigs[generation][0];
        config_pop->population_size = populationSize_;
        config_pop->mating_distribution.push_back(MatingDistribution::Entry(1, 0,0));
    }

    // write/read test 

    ostringstream oss;
    oss << populationConfigs;

    vector<Population::Configs> populationConfigs2;
    unit_assert(populationConfigs != populationConfigs2);

    istringstream iss(oss.str());
    iss >> populationConfigs2;

    // check automatic id_offset setting
    
    unit_assert(populationConfigs2.size() == generation_count);
    unit_assert(populationConfigs2[0].size() == 3);

    unit_assert(populationConfigs2[0][0].id_offset == 0);
    unit_assert(populationConfigs2[0][1].id_offset == 0);
    unit_assert(populationConfigs2[0][2].id_offset == populationSize_*2);

    populationConfigs2[0][2].id_offset = 0; // fix up for equality comparison
    unit_assert(populationConfigs == populationConfigs2);
}


void test_Population_IO()
{
    if (os_) *os_ << "test_Population_IO()\n";

    Population::Config config;
    config.population_size = 10;
    config.chromosome_pair_count = 1;
    Population_Organisms p;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    p.create_organisms(config, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    unit_assert(p.organisms().size() == 10);
    if (os_) *os_ << "p:\n" << p << endl;

    ostringstream oss;
    oss << p;

    // test reading into Population_Organisms

    Population_Organisms q;
    unit_assert(p != q);

    istringstream iss(oss.str());
    iss >> q;
    if (os_) *os_ << "q:\n" << q << endl;
    unit_assert(p == q);

    if (os_) *os_ << endl;

    // test reading into Population_ChromosomePairs

    Population_ChromosomePairs r;
    unit_assert(p != r);

    iss.seekg(0);
    iss >> r;
    if (os_) *os_ << "r:\n" << r << endl;
    unit_assert(p == r);

    if (os_) *os_ << endl;
}


void test_Population_IO_Binary()
{
    if (os_) *os_ << "test_Population_IO_Binary()\n";

    Population::Config config;
    config.population_size = 10;
    config.chromosome_pair_count = 1;
    Population_Organisms p;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    p.create_organisms(config, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    unit_assert(p.organisms().size() == 10);
    if (os_) *os_ << "p:\n" << p << endl;

    ostringstream oss;
    p.write_binary(oss);

    // test reading into Population_Organisms

    Population_Organisms q;
    unit_assert(p != q);

    istringstream iss(oss.str());
    q.read_binary(iss);
    if (os_) *os_ << "q:\n" << q << endl;
    unit_assert(p == q);

    // test reading into Population_ChromosomePairs

    Population_ChromosomePairs r;
    unit_assert(p != r);

    iss.seekg(0);
    r.read_binary(iss);
    if (os_) *os_ << "r:\n" << r << endl;
    unit_assert(p == r);

    if (os_) *os_ << endl;
}


void demo_Population_mutate()
{
    if (os_) *os_ << "demo_Population_mutate()\n";

    Population::Config config;
    config.population_size = 10;
    config.chromosome_pair_count = 1;
    //Population_Organisms p;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    Population_ChromosomePairs p;
    p.create_organisms(config, PopulationPtrs(), PopulationDataPtrs(), rpgs);
    unit_assert(p.population_size() == 10);
    if (os_) *os_ << "p:\n" << p << endl;

    size_t individual_index = 6;
    size_t chromosome_pair_index = 0;
    unsigned int position = 1000;
    unsigned int new_id = 666;

    HaplotypeChunk& chunk = *(p.chromosome_pair_range(individual_index).begin()+chromosome_pair_index)
        ->first.find_haplotype_chunk(position);

    chunk.id = new_id;

    if (os_) *os_ << "mutation\n\np:\n" << p << endl;
}


void test_Population_create()
{
    Population::Configs configs_gen0(2);
    configs_gen0[0].population_size = 10;
    configs_gen0[0].chromosome_pair_count = 1;
    configs_gen0[0].id_offset = 0;
    configs_gen0[1].population_size = 10;
    configs_gen0[1].chromosome_pair_count = 1;
    configs_gen0[1].id_offset = 100;

    PopulationPtrs dummy;
    PopulationDataPtrs dummy_datas;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    PopulationPtrsPtr populations_gen0 = Population::create_populations(configs_gen0,
                                                                        dummy,
                                                                        dummy_datas,
                                                                        rpgs);
    if (os_)
    {
        for (PopulationPtrs::const_iterator it=populations_gen0->begin();
            it!=populations_gen0->end(); ++it)
            *os_ << **it << endl;
    }

    const Population& pop0 = *populations_gen0->at(0);
    const Population& pop1 = *populations_gen0->at(1);

    // validate id assignment:
    //   pop0: 0, 1, ..., 19
    //   pop1: 100, 101, ..., 119

    for (size_t i=0; i<10; ++i)
    {
        // check ith individual in pop0
        const ChromosomePair& p = *pop0.chromosome_pair_range(i).begin();
        unit_assert(p.first.haplotype_chunks().size() == 1);
        unit_assert(p.first.haplotype_chunks().front().id == 2*i);
        unit_assert(p.second.haplotype_chunks().size() == 1);
        unit_assert(p.second.haplotype_chunks().front().id == 2*i + 1);

        // check ith individual in pop1
        const ChromosomePair& p2 = *pop1.chromosome_pair_range(i).begin();
        unit_assert(p2.first.haplotype_chunks().size() == 1);
        unit_assert(p2.first.haplotype_chunks().front().id == 100 + 2*i);
        unit_assert(p2.second.haplotype_chunks().size() == 1);
        unit_assert(p2.second.haplotype_chunks().front().id == 100 + 2*i + 1);
    }

    PopulationDataPtrs population_datas;
    population_datas.push_back(PopulationDataPtr(new PopulationData));
    population_datas.push_back(PopulationDataPtr(new PopulationData));
    population_datas[0]->population_size = 10;
    population_datas[1]->population_size = 10;

    double fitness_raw[] = {1,1,1,1,1,0,0.0,0,0};
    double fitness_special_raw[] = {1,0,0,0,0,0,0.0,0,0};

    DataVectorPtr fitness(new DataVector(10));
    copy(fitness_raw, fitness_raw + sizeof(fitness_raw)/sizeof(double), fitness->begin());
    DataVectorPtr fitness_special(new DataVector(10));
    copy(fitness_special_raw, fitness_special_raw + sizeof(fitness_special_raw)/sizeof(double), fitness_special->begin());

    (*population_datas[0]->trait_values)["fitness"] = fitness;
    (*population_datas[0]->trait_values)["fitness_special"] = fitness_special;
    (*population_datas[1]->trait_values)["fitness"] = fitness;
    (*population_datas[1]->trait_values)["fitness_special"] = fitness_special;

    Population::Configs configs_gen1(1);
    Population::Config& config = configs_gen1[0];
    config.population_size = 10;
    config.chromosome_pair_count = 1;

    // select mom from first 5 individuals in pop0, dad from first 5 in pop1

    config.mating_distribution.default_fitness_function = "fitness";
    config.mating_distribution.push_back(MatingDistribution::Entry(1.0, 0, 1));

    PopulationPtrsPtr populations_gen1 = Population::create_populations(configs_gen1,
                                                                        *populations_gen0,
                                                                        population_datas,
                                                                        rpgs);
    const Population& pop_gen1 = *populations_gen1->at(0);
    if (os_) *os_ << "gen1:\n" << pop_gen1 << endl;

    unit_assert(pop_gen1.population_size() == 10);
    for (size_t i=0; i<10; ++i)
    {
        const ChromosomePair& p = *pop_gen1.chromosome_pair_range(i).begin();
        unit_assert(p.first.haplotype_chunks().size() == 1);
        unit_assert(p.first.haplotype_chunks().front().id < 10);
        unit_assert(p.second.haplotype_chunks().size() == 1);
        unit_assert(p.second.haplotype_chunks().front().id > 20 &&
                    p.second.haplotype_chunks().front().id < 110);
    }

    // fitness_special: mom is first individual in pop1, dad in first five in pop0

    config.mating_distribution = MatingDistribution(); // clear
    config.mating_distribution.default_fitness_function = "fitness";
    config.mating_distribution.push_back(
        MatingDistribution::Entry(1.0, 1, "fitness_special", 0, ""));

    populations_gen1 = Population::create_populations(configs_gen1,
                                                      *populations_gen0,
                                                      population_datas,
                                                      rpgs);

    const Population& pop_gen1_special = *populations_gen1->at(0);
    if (os_) *os_ << "gen1 special:\n" << pop_gen1_special << endl;

    unit_assert(pop_gen1_special.population_size() == 10);

    for (size_t i=0; i<10; ++i)
    {
        const ChromosomePair& p = *pop_gen1_special.chromosome_pair_range(i).begin();
        unit_assert(p.first.haplotype_chunks().size() == 1);
        unit_assert(p.first.haplotype_chunks().front().id > 20 &&
                    (p.first.haplotype_chunks().front().id == 100 ||
                    p.first.haplotype_chunks().front().id == 101));
        unit_assert(p.second.haplotype_chunks().size() == 1);
        unit_assert(p.second.haplotype_chunks().front().id < 10);
    }
}


void test()
{
    test_MatingDistribution();
    test_MatingDistribution_validate_entries();
    test_Population_Config_IO();
    test_generations_IO();
    test_Population_IO();
    test_Population_IO_Binary();
    demo_Population_mutate();
    test_Population_create();
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


