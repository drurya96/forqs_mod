//
// MutationGeneratorImplementationTest.cpp
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


#include "MutationGeneratorImplementation.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "Population_ChromosomePairs.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>
#include <numeric>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_MutationGenerator_SingleLocus()
{
    if (os_) *os_ << "test_MutationGenerator_SingleLocus()\n\n";

    LocusPtr locus1(new Locus("locus1", 0, 500));

    Configurable::Registry registry;
    registry["locus1"] = locus1;

    Parameters parameters;
    parameters.insert_name_value("locus", "locus1");
    parameters.insert_name_value("mu", .5);

    MutationGenerator_SingleLocus mg("mg");
    mg.configure(parameters, registry);

    Parameters parameters_out = mg.parameters();

    if (os_)
    {
        *os_ << "parameters:\n\n" << parameters << "\n";
        *os_ << "parameters_out:\n\n" << parameters_out << "\n";
        *os_ << "configuration:\n\n";
        set<string> ids_written;
        mg.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    Population::Config config;
    config.population_size = 100;
    config.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    Population_ChromosomePairs population;
    population.create_organisms(config, PopulationPtrs(), PopulationDataPtrs(), rpgs);

    if (os_)
        *os_ << "population size: " << population.population_size() << endl;

    double sum_mutant_count = 0;
    double sum_which = 0;
    double sum_even = 0;
    const size_t generation_count = 100;
    const size_t population_index = 0;

    for (size_t generation_index=0; generation_index<generation_count; ++generation_index)
    {
        MutationGenerator::MutationInfos mutation_infos = 
            mg.generate_mutations(population, generation_index, population_index);

        sum_mutant_count += mutation_infos.size();

        for (MutationGenerator::MutationInfos::const_iterator it=mutation_infos.begin();
             it!=mutation_infos.end(); ++it)
        {
            sum_which += it->which;
            sum_even += (it->individual_index % 2 == 0);
        }

        /*
        if (os_)
        {
            copy(mutation_infos.begin(), mutation_infos.end(), 
                 ostream_iterator<MutationGenerator::MutationInfo>(*os_, "\n"));
            *os_ << endl;
        }
        */
    }

    double average_mutant_count = sum_mutant_count / generation_count;
    double average_which = sum_which / sum_mutant_count;
    double average_even = sum_even / sum_mutant_count;

    if (os_)
        *os_ << "average_mutant_count: " << average_mutant_count << endl
             << "average_which: " << average_which << endl
             << "average_even: " << average_even << endl;

    unit_assert_equal(average_mutant_count, config.population_size, .1); // 2 * population_size chromosomes, mu == .5
    unit_assert_equal(average_which, .5, .005);
    unit_assert_equal(average_even, .5, .005);
}


void test_MutationGenerator_Regions()
{
    if (os_) *os_ << "test_MutationGenerator_Regions()\n\n";

    LocusPtr locus1(new Locus("locus1", 0, 100000));
    LocusPtr locus2(new Locus("locus2", 1, 200000));
    LocusPtr locus3(new Locus("locus3", 2, 300000));

    TrajectoryPtr mu1(new Trajectory_Constant("mu1", 1e-4));
    TrajectoryPtr mu2(new Trajectory_Constant("mu2", 1e-5));
    TrajectoryPtr mu3(new Trajectory_Constant("mu3", 1e-6));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;
    registry["locus3"] = locus3;
    registry["mu1"] = mu1;
    registry["mu2"] = mu2;
    registry["mu3"] = mu3;

    Parameters parameters;
    parameters.insert_name_value("locus:length:rate", "locus1 100 mu1");
    parameters.insert_name_value("locus:length:rate", "locus2 1000 mu2");
    parameters.insert_name_value("locus:length:rate", "locus3 10000 mu3");

    MutationGenerator_Regions mg("id_dummy");
    mg.configure(parameters, registry);

    Parameters parameters_out = mg.parameters();

    if (os_)
    {
        *os_ << "parameters:\n\n" << parameters << "\n";
        *os_ << "parameters_out:\n\n" << parameters_out << "\n";
        *os_ << "configuration:\n\n";
        set<string> ids_written;
        mg.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    Population::Config config;
    config.population_size = 1000;
    config.chromosome_pair_count = 1;

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_Trivial("rpg")));
    rpgs.push_back(rpgs.front());

    Population_ChromosomePairs population;
    population.create_organisms(config, PopulationPtrs(), PopulationDataPtrs(), rpgs);

    // locus1:  100 sites * 1e-4 mutations/site/chromosome/generation * (2*1000) chromosomes == 20 mutations/generation 
    // locus2:  1000 sites * 1e-5 mutations/site/chromosome/generation * (2*1000) chromosomes == 20 mutations/generation 
    // locus3:  10000 sites * 1e-6 mutations/site/chromosome/generation * (2*1000) chromosomes == 20 mutations/generation 

    if (os_)
        *os_ << "population size: " << population.population_size() << endl;

    vector<double> chromosome_sums(3);
    double sum_mutant_count = 0;
    double sum_which = 0;
    const size_t generation_count = 30;
    const size_t population_index = 0;

    for (size_t generation_index=0; generation_index<generation_count; ++generation_index)
    {
        MutationGenerator::MutationInfos mutation_infos = 
            mg.generate_mutations(population, generation_index, population_index);

        sum_mutant_count += mutation_infos.size();

        for (MutationGenerator::MutationInfos::const_iterator it=mutation_infos.begin();
             it!=mutation_infos.end(); ++it)
        {
            ++chromosome_sums[it->locus.chromosome_pair_index];
            sum_which += it->which;
        }

        /*
        if (os_)
        {
            copy(mutation_infos.begin(), mutation_infos.end(), 
                 ostream_iterator<MutationGenerator::MutationInfo>(*os_, "\n"));
            *os_ << endl;
        }
        */
    }

    unit_assert(sum_mutant_count == accumulate(chromosome_sums.begin(), chromosome_sums.end(), 0.));

    vector<double> average_mutant_counts(3);
    transform(chromosome_sums.begin(), chromosome_sums.end(), average_mutant_counts.begin(),
              bind2nd(divides<double>(), generation_count));

    double average_which = sum_which / sum_mutant_count;

    if (os_)
    {
        *os_ << "average_mutant_counts: ";
        copy(average_mutant_counts.begin(), average_mutant_counts.end(), ostream_iterator<double>(*os_, " "));
        *os_ << endl;
        *os_ << "average_which: " << average_which << endl;
    }

    for (vector<double>::const_iterator it=average_mutant_counts.begin(); it!=average_mutant_counts.end(); ++it)
        unit_assert_equal(*it, 20, 3);
    unit_assert_equal(average_which, .5, 2e-2);
}


void test()
{
    test_MutationGenerator_SingleLocus();
    test_MutationGenerator_Regions();
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


