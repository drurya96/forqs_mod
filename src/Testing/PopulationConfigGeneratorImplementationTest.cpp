//
// PopulationConfigGeneratorImplementationTest.cpp
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


#include "PopulationConfigGeneratorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_Configurable_PopulationConfigGenerator_File()
{
    if (os_) *os_ << "test_Configurable_PopulationConfigGenerator_File()\n";

    Parameters parameters_in;
    parameters_in.insert_name_value("filename", "../examples/popconfig_10k.txt");

    PopulationConfigGenerator_File pop_config_generator("dummy_id");
    unit_assert(pop_config_generator.generation_count() == 0);

    Configurable::Registry registry;
    pop_config_generator.configure(parameters_in, registry);

    unit_assert(pop_config_generator.generation_count() == 7);
    Population::Configs configs = pop_config_generator.population_configs(0, PopulationDataPtrs());
    unit_assert(configs.size() == 3);
    unit_assert(configs[0].population_size == 0);
    unit_assert(configs[1].population_size == 10000);
    unit_assert(configs[2].population_size == 10000);

    for (size_t i=1; i<8; ++i)
    {
        configs = pop_config_generator.population_configs(i, PopulationDataPtrs());
        unit_assert(configs.size() == 1);
        unit_assert(configs[0].population_size == 10000);
    }

    Parameters parameters_out = pop_config_generator.parameters();
    unit_assert(parameters_in == parameters_out);

    if (os_) *os_ << "min_unused_id: " << pop_config_generator.min_unused_id() << endl;
    unit_assert(pop_config_generator.min_unused_id() == 40000);
}


void test_Configurable_PopulationConfigGenerator_ConstantSize()
{
    if (os_) *os_ << "test_Configurable_PopulationConfigGenerator_ConstantSize()\n";

    Parameters parameters_in;
    parameters_in.insert_name_value("chromosome_pair_count", 12);
    parameters_in.insert_name_value("population_size", 34);
    parameters_in.insert_name_value("population_count", 56);
    parameters_in.insert_name_value("generation_count", 78);
    parameters_in.insert_name_value("id_offset_step", 0);
    parameters_in.insert_name_value("chromosome_lengths", "1 2 3 4 5 6 7 8 9 10 11 12 ");
    
    if (os_) *os_ << "parameters_in:\n" << parameters_in << endl;;

    PopulationConfigGenerator_ConstantSize pop_config_generator("dummy_id");

    Configurable::Registry registry;
    pop_config_generator.configure(parameters_in, registry);
    Parameters parameters_out = pop_config_generator.parameters();

    if (os_) *os_ << "parameters_out:\n" << parameters_out << endl;;

    unit_assert(parameters_in == parameters_out);

    if (os_) *os_ << "min_unused_id: " << pop_config_generator.min_unused_id() << endl;
    unit_assert(pop_config_generator.min_unused_id() == 2*34*56);

    unit_assert(pop_config_generator.chromosome_lengths().size() == 12);
}


void test_PopulationConfigGenerator_LinearSteppingStone()
{
    if (os_) *os_ << "test_PopulationConfigGenerator_LinearSteppingStone()\n";

    string id_population_size_trajectory("id_population_size_trajectory");
    string id_migration_rate_trajectory_default("id_migration_rate_trajectory_default");
    string id_migration_rate_trajectory_123("id_123");

    TrajectoryPtr population_size_trajectory(new Trajectory_Linear(id_population_size_trajectory, 2000, 1000));
    TrajectoryPtr migration_rate_trajectory_default(new Trajectory_Constant(id_migration_rate_trajectory_default, .05));
    TrajectoryPtr migration_rate_trajectory_123(new Trajectory_Constant(id_migration_rate_trajectory_123, .123));

    Configurable::Registry registry;
    registry[id_population_size_trajectory] = population_size_trajectory;
    registry[id_migration_rate_trajectory_default] = migration_rate_trajectory_default;
    registry[id_migration_rate_trajectory_123] = migration_rate_trajectory_123;

    const size_t population_count = 3;
    const size_t generation_count = 10;
    const size_t id_offset_step = 10000;

    Parameters parameters;
    parameters.insert_name_value("chromosome_pair_count", 1); 
    parameters.insert_name_value("population_count", population_count);
    parameters.insert_name_value("generation_count", generation_count);
    parameters.insert_name_value("id_offset_step", id_offset_step);
    parameters.insert_name_value("population_size", id_population_size_trajectory);
    parameters.insert_name_value("migration_rate_default", id_migration_rate_trajectory_default);
    parameters.insert_name_value("migration_rate:from:to", "id_123 2 3"); // 1-based
    parameters.insert_name_value("migration_rate:from:to", "id_123 1 2"); // 1-based

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_LinearSteppingStone("id_pcg"));
    pcg->configure(parameters, registry);

    Parameters parameters_out = pcg->parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
        *os_ << "configuration:\n";
        set<string> written;
        pcg->write_configuration(*os_, written);
    }

    unit_assert(parameters == parameters_out);


    // test output

    if (os_) *os_ << "min_unused_id: " << pcg->min_unused_id() << endl;
    unit_assert(pcg->min_unused_id() == (population_count-1) * id_offset_step + 2*1000); // last pop size == 1000

    vector<MatingDistribution> mds(3);

    mds[0].push_back(MatingDistribution::Entry(.95, 0, 0));
    mds[0].push_back(MatingDistribution::Entry(.05, 1, 1));

    mds[1].push_back(MatingDistribution::Entry(.827, 1, 1));
    mds[1].push_back(MatingDistribution::Entry(.123, 0, 0));
    mds[1].push_back(MatingDistribution::Entry(.05, 2, 2));

    mds[2].push_back(MatingDistribution::Entry(.877, 2, 2));
    mds[2].push_back(MatingDistribution::Entry(.123, 1, 1));

    if (os_) *os_ << "popconfigs:\n";

    for (size_t generation_index=0; generation_index<=generation_count; ++generation_index)
    {
        Population::Configs popconfigs = pcg->population_configs(generation_index, PopulationDataPtrs());

        if (os_) 
        {
            *os_ << "generation " << generation_index << endl;
            copy(popconfigs.begin(), popconfigs.end(), ostream_iterator<Population::Config>(*os_, "\n"));
            *os_ << endl;
        }

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            unit_assert(popconfigs[population_index].population_size == 1000 + 2000*generation_index);
            if (generation_index > 0)
            {
                unit_assert(popconfigs[population_index].mating_distribution == mds[population_index]);
            }
            else
            {
                unit_assert(popconfigs[population_index].mating_distribution.empty());
                unit_assert(popconfigs[population_index].id_offset == id_offset_step * population_index);
            }
        }
    }
}


void test_PopulationConfigGenerator_Island()
{
    if (os_) *os_ << "test_PopulationConfigGenerator_Island()\n";

    string id_population_size_trajectory("id_population_size_trajectory");
    string id_migration_rate_trajectory_default("id_migration_rate_trajectory_default");
    string id_migration_rate_trajectory_123("id_123");

    TrajectoryPtr population_size_trajectory(new Trajectory_Linear(id_population_size_trajectory, 2000, 1000));
    TrajectoryPtr migration_rate_trajectory_default(new Trajectory_Constant(id_migration_rate_trajectory_default, .05));
    TrajectoryPtr migration_rate_trajectory_123(new Trajectory_Constant(id_migration_rate_trajectory_123, .123));

    Configurable::Registry registry;
    registry[id_population_size_trajectory] = population_size_trajectory;
    registry[id_migration_rate_trajectory_default] = migration_rate_trajectory_default;
    registry[id_migration_rate_trajectory_123] = migration_rate_trajectory_123;

    const size_t population_count = 3;
    const size_t generation_count = 10;
    const size_t id_offset_step = 10000;

    Parameters parameters;
    parameters.insert_name_value("chromosome_pair_count", 1); 
    parameters.insert_name_value("population_count", population_count);
    parameters.insert_name_value("generation_count", generation_count);
    parameters.insert_name_value("id_offset_step", id_offset_step);
    parameters.insert_name_value("population_size", id_population_size_trajectory);
    parameters.insert_name_value("migration_rate_default", id_migration_rate_trajectory_default);
    parameters.insert_name_value("migration_rate:from:to", "id_123 2 3"); // 1-based
    parameters.insert_name_value("migration_rate:from:to", "id_123 1 2"); // 1-based

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_Island("id_pcg"));
    pcg->configure(parameters, registry);

    Parameters parameters_out = pcg->parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
        *os_ << "configuration:\n";
        set<string> written;
        pcg->write_configuration(*os_, written);
    }

    unit_assert(parameters == parameters_out);


    // test output

    if (os_) *os_ << "min_unused_id: " << pcg->min_unused_id() << endl;
    unit_assert(pcg->min_unused_id() == (population_count-1) * id_offset_step + 2*1000); // last pop size == 1000

    vector<MatingDistribution> mds(3);

    mds[0].push_back(MatingDistribution::Entry(.90, 0, 0));
    mds[0].push_back(MatingDistribution::Entry(.05, 1, 1));
    mds[0].push_back(MatingDistribution::Entry(.05, 2, 2));

    mds[1].push_back(MatingDistribution::Entry(.827, 1, 1));
    mds[1].push_back(MatingDistribution::Entry(.123, 0, 0));
    mds[1].push_back(MatingDistribution::Entry(.05, 2, 2));

    mds[2].push_back(MatingDistribution::Entry(.827, 2, 2));
    mds[2].push_back(MatingDistribution::Entry(.05, 0, 0));
    mds[2].push_back(MatingDistribution::Entry(.123, 1, 1));

    if (os_) *os_ << "popconfigs:\n";

    for (size_t generation_index=0; generation_index<=generation_count; ++generation_index)
    {
        Population::Configs popconfigs = pcg->population_configs(generation_index, PopulationDataPtrs());

        if (os_) 
        {
            *os_ << "generation " << generation_index << endl;
            copy(popconfigs.begin(), popconfigs.end(), ostream_iterator<Population::Config>(*os_, "\n"));
            *os_ << endl;
        }

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            unit_assert(popconfigs[population_index].population_size == 1000 + 2000*generation_index);
            if (generation_index > 0)
            {
                unit_assert(popconfigs[population_index].mating_distribution == mds[population_index]);
            }
            else
            {
                unit_assert(popconfigs[population_index].mating_distribution.empty());
                unit_assert(popconfigs[population_index].id_offset == id_offset_step * population_index);
            }
        }
    }
}


void test()
{
    test_Configurable_PopulationConfigGenerator_File();
    test_Configurable_PopulationConfigGenerator_ConstantSize();
    test_PopulationConfigGenerator_LinearSteppingStone();
    test_PopulationConfigGenerator_Island();
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


