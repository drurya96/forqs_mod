//
// PopulationConfigGenerator.cpp
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


#include "PopulationConfigGenerator.hpp"
#include "boost/lambda/casts.hpp"
#include <fstream>
#include <stdexcept>
#include <numeric>


using namespace std;
using namespace boost::lambda; // for _1 and ll_static_cast


//
// PopulationConfigGenerator
//


const vector<string> default_fitness_functions_(1); // single empty string


PopulationConfigGenerator::PopulationConfigGenerator(const string& id)
:   Configurable(id),
    generation_count_(0),
    population_count_(0),
    id_offset_step_(0),
    chromosome_pair_count_(0),
    fitness_functions_(default_fitness_functions_)
{}


unsigned int PopulationConfigGenerator::min_unused_id() const
{
    Population::Configs configs = population_configs(0, PopulationDataPtrs());
    unsigned int result = 0;
    for (Population::Configs::const_iterator it=configs.begin(); it!=configs.end(); ++it)
        if (result < it->id_offset + 2*it->population_size)
            result = it->id_offset + 2*it->population_size;
    return result;
}


string PopulationConfigGenerator::class_name() const
{
    cerr << "[PopulationConfigGenerator] Warning: virtual class_name() has not been defined in derived class.\n";
    return "PopulationConfigGenerator";
}


Parameters PopulationConfigGenerator::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("generation_count", generation_count_);
    parameters.insert_name_value("population_count", population_count_);
    parameters.insert_name_value("id_offset_step", id_offset_step_);
    parameters.insert_name_value("chromosome_pair_count", chromosome_pair_count_);
    if (fitness_functions_ != default_fitness_functions_)
        parameters.insert_name_value_vector("fitness_function", fitness_functions_);
    if (!chromosome_lengths_.empty())
        parameters.insert_name_value_vector("chromosome_lengths", chromosome_lengths_);
    return parameters;
}


void PopulationConfigGenerator::configure(const Parameters& parameters, const Registry& registry)
{
    generation_count_ = parameters.value<size_t>("generation_count", 0);
    population_count_ = parameters.value<size_t>("population_count", 1);
    id_offset_step_ = parameters.value<size_t>("id_offset_step", 0);
    chromosome_pair_count_ = parameters.value<size_t>("chromosome_pair_count", 1);

    fitness_functions_ = parameters.value_vector<string>("fitness_function", default_fitness_functions_);

    if (parameters.count("chromosome_lengths"))
    {
        vector<double> lengths = parameters.value_vector<double>("chromosome_lengths");
        if (lengths.size() != chromosome_pair_count_)
            throw runtime_error("[PopulationConfigGenerator] Chromosome length count does not match chromosome pair count.");
        chromosome_lengths_.clear();
        transform(lengths.begin(), lengths.end(), back_inserter(chromosome_lengths_), ll_static_cast<unsigned int>(_1));
    }
}


