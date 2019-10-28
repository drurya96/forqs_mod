//
// PopulationConfigGeneratorImplementation.cpp
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
#include "boost/lambda/casts.hpp"
#include <fstream>
#include <stdexcept>
#include <numeric>


using namespace std;
using namespace boost::lambda; // for _1 and ll_static_cast


//
// PopulationConfigGenerator_File
//


PopulationConfigGenerator_File::PopulationConfigGenerator_File(const string& id, const string& filename)
:   PopulationConfigGenerator(id), filename_(filename)
{
    if (!filename_.empty()) read_file();
}


Population::Configs PopulationConfigGenerator_File::population_configs(size_t generation_index,
    const PopulationDataPtrs& population_datas) const
{
    if (generation_index >= population_configs_.size())
        throw runtime_error("[PopulationConfigGenerator_File::population_configs()] Invalid generation index.");

    return population_configs_[generation_index];
}


Parameters PopulationConfigGenerator_File::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("filename", filename_);
    if (!chromosome_lengths_.empty())
        parameters.insert_name_value_vector("chromosome_lengths", chromosome_lengths_);
    return parameters;
}


void PopulationConfigGenerator_File::configure(const Parameters& parameters, const Registry& registry)
{
    filename_ = parameters.value<string>("filename");
    read_file();

    if (parameters.count("chromosome_lengths"))
    {
        vector<double> lengths = parameters.value_vector<double>("chromosome_lengths");
        if (lengths.size() != chromosome_pair_count_)
            throw runtime_error("[PopulationConfigGenerator] Chromosome length count does not match chromosome pair count.");
        chromosome_lengths_.clear();
        transform(lengths.begin(), lengths.end(), back_inserter(chromosome_lengths_), ll_static_cast<unsigned int>(_1));
    }
}


void PopulationConfigGenerator_File::read_file()
{
    population_configs_.clear();

    ifstream is(filename_.c_str());
    if (!is)
        throw runtime_error(("[PopulationConfigGenerator_File] Unable to open file " + filename_).c_str());
    
    is >> population_configs_;

    // set PCG protected members

    generation_count_ = population_configs_.empty() ? 0 : population_configs_.size() - 1;

    if (!population_configs_.empty())
    {
        const Population::Configs& popconfigs_gen0 = population_configs_.front();
        population_count_ = popconfigs_gen0.size();

        if (!popconfigs_gen0.empty())
            chromosome_pair_count_ = popconfigs_gen0.front().chromosome_pair_count;
    }
}


//
// PopulationConfigGenerator_ConstantSize
//


Population::Configs PopulationConfigGenerator_ConstantSize::population_configs(
    size_t generation_index, const PopulationDataPtrs& population_datas) const
{
    if (generation_index > generation_count_)
        throw runtime_error("[PopulationConfigGenerator_ConstantSize] Invalid generation index.");

    Population::Configs popconfigs(population_count_);
    for (size_t i=0; i<population_count_; ++i)
    {
        Population::Config& popconfig = popconfigs[i];

        popconfig.population_size = population_size_;
        popconfig.chromosome_pair_count = chromosome_pair_count_;

        if (generation_index == 0)
        {
            unsigned int step = id_offset_step_ > 0 ? id_offset_step_ : 2 * population_size_;            
            popconfig.id_offset = i * step;
        }

        if (generation_index > 0)
            popconfig.mating_distribution.push_back(MatingDistribution::Entry(1, i, i));

        if (fitness_functions_.size() == 1)
            popconfig.mating_distribution.default_fitness_function = fitness_functions_.front();
        else if (i < fitness_functions_.size())
            popconfig.mating_distribution.default_fitness_function = fitness_functions_[i];
        else
            throw runtime_error("[PopulationConfigGenerator_ConstantSize] Invalid fitness function specification.");
    }
    
    return popconfigs;
}


Parameters PopulationConfigGenerator_ConstantSize::parameters() const
{
    Parameters parameters = PopulationConfigGenerator::parameters();
    parameters.insert_name_value("population_size", population_size_);
    return parameters;
}


void PopulationConfigGenerator_ConstantSize::configure(const Parameters& parameters, const Registry& registry)
{
    PopulationConfigGenerator::configure(parameters, registry);
    population_size_ = parameters.value<size_t>("population_size");

    if (id_offset_step_ != 0 && id_offset_step_ < 2 * population_size_)
        throw runtime_error("[PopulationConfigGenerator_ConstantSize] id_offset_step < 2*population_size");
}


//
// MigrationRateTrajectoryInfo
//


MigrationRateTrajectoryInfo::MigrationRateTrajectoryInfo(
    const string& configuration, 
    const Configurable::Registry& registry)
{
    string id_trajectory;
    istringstream iss(configuration);

    size_t population_from = 0, population_to = 0;
    iss >> id_trajectory >> population_from >> population_to; // 1-based
        
    if (population_from == 0 || population_to == 0)
        throw runtime_error("[MigrationRateTrajectoryInfo] Invalid population == 0.");

    population_index_from = population_from - 1;
    population_index_to = population_to - 1;

    if (id_trajectory.empty())
    {
        ostringstream oss;
        oss << "[PopulationConfigGenerator_LinearSteppingStone] Required format:\n"
            << "\tmigration_rate_default <id_trajectory>\n"
            << "or:\n"
            << "\tmigration_rate:from:to <id_trajectory> <int_population_from> <int_population_to>\n";
        throw runtime_error(oss.str().c_str());
    }

    trajectory = registry.get<Trajectory>(id_trajectory);    

    if (!trajectory.get())
        throw runtime_error("[PopulationConfigGenerator_LinearSteppingStone] Null trajectory.");
}


string MigrationRateTrajectoryInfo::configuration() const
{
    ostringstream oss;
    oss << trajectory->object_id() << " "
        << population_index_from + 1 << " "
        << population_index_to + 1;
    return oss.str();
}


//
// PopulationConfigGenerator_IslandBase
//


PopulationConfigGenerator_IslandBase::PopulationConfigGenerator_IslandBase(const string& id)
:   PopulationConfigGenerator(id)
{}


Population::Configs PopulationConfigGenerator_IslandBase::population_configs(
    size_t generation_index, const PopulationDataPtrs& population_datas) const
{
    if (generation_index > generation_count_)
        throw runtime_error("[PopulationConfigGenerator_IslandBase] Invalid generation index.");

    Population::Configs popconfigs(population_count_);

    unsigned int current_id_offset = 0;

    for (size_t population_index=0; population_index<population_count_; ++population_index)
    {
        Population::Config& popconfig = popconfigs[population_index];

        popconfig.chromosome_pair_count = chromosome_pair_count_;
        popconfig.population_size = size_t(round(population_size_trajectory_->value(generation_index, population_index)));

        if (fitness_functions_.size() == 1)
            popconfig.mating_distribution.default_fitness_function = fitness_functions_.front();
        else if (population_index < fitness_functions_.size())
            popconfig.mating_distribution.default_fitness_function = fitness_functions_[population_index];
        else
            throw runtime_error("[PopulationConfigGenerator_ConstantSize] Invalid fitness function specification.");

        if (generation_index == 0)
        {
            if (id_offset_step_ > 0)
            {
                if (id_offset_step_ < 2 * popconfig.population_size)
                    throw runtime_error("[PopulationConfigGenerator_IslandBase] Error: id_offset_step < 2 * initial population size.");

                popconfig.id_offset = population_index * id_offset_step_;
            }
            else
            {
                popconfig.id_offset = current_id_offset;
                current_id_offset += 2 * popconfig.population_size;
            }
        }
        else // if (generation_index > 0)
        {
            vector<double> immigration_rates(population_count_);

            for (size_t from_index=0; from_index<population_count_; ++from_index)
            {
                immigration_rates[from_index] = (from_index == population_index) ? 0 :
                    migration_rate_trajectories_[population_index][from_index]->value(generation_index, population_index);
            }

            double total_immigration_rate = accumulate(immigration_rates.begin(), immigration_rates.end(), 0.);

            if (total_immigration_rate > 1.0)
            {
                ostringstream oss;
                oss << "[PopulationConfigGenerator_IslandBase] Error: total immigration rate "
                    << total_immigration_rate << " > 1.0 "
                    << "in generation " << generation_index << " population " << population_index + 1;
                throw runtime_error(oss.str().c_str());
            }

            popconfig.mating_distribution.push_back(
                 MatingDistribution::Entry(1 - total_immigration_rate, 
                                           population_index, 
                                           population_index));
    
            for (size_t from_index=0; from_index<population_count_; ++from_index)
                if (immigration_rates[from_index] > 0.) 
                    popconfig.mating_distribution.push_back(
                        MatingDistribution::Entry(immigration_rates[from_index], 
                                                  from_index, 
                                                  from_index));
        }
    }
    
    return popconfigs;
}


Parameters PopulationConfigGenerator_IslandBase::parameters() const
{
    Parameters parameters = PopulationConfigGenerator::parameters();

    if (population_size_trajectory_.get())
        parameters.insert_name_value("population_size", population_size_trajectory_->object_id());
    if (migration_rate_trajectory_default_.get())
        parameters.insert_name_value("migration_rate_default", migration_rate_trajectory_default_->object_id());

    for (MigrationRateTrajectoryInfos::const_iterator info=migration_rate_trajectory_infos_.begin();
        info!=migration_rate_trajectory_infos_.end(); ++info)
        parameters.insert_name_value("migration_rate:from:to", info->configuration());

    return parameters;
}


void PopulationConfigGenerator_IslandBase::configure(const Parameters& parameters, const Registry& registry)
{
    PopulationConfigGenerator::configure(parameters, registry);

    population_size_trajectory_ = 
        registry.get<Trajectory>(parameters.value<string>("population_size"));

    if (parameters.count("migration_rate_default"))
        migration_rate_trajectory_default_ = 
            registry.get<Trajectory>(parameters.value<string>("migration_rate_default"));
    
    vector<string> configurations = parameters.values<string>("migration_rate:from:to");

    for (vector<string>::const_iterator configuration=configurations.begin(); configuration!=configurations.end(); ++configuration)
        migration_rate_trajectory_infos_.push_back(MigrationRateTrajectoryInfo(*configuration, registry));

    initialize_trajectories();
}


void PopulationConfigGenerator_IslandBase::write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const
{
    if (population_size_trajectory_.get())
        population_size_trajectory_->write_configuration(os, ids_written);

    if (migration_rate_trajectory_default_.get())
        migration_rate_trajectory_default_->write_configuration(os, ids_written);

    for (MigrationRateTrajectoryInfos::const_iterator info=migration_rate_trajectory_infos_.begin();
        info!=migration_rate_trajectory_infos_.end(); ++info)
        info->trajectory->write_configuration(os, ids_written);
}


void PopulationConfigGenerator_IslandBase::initialize_trajectories()
{
    if (!population_size_trajectory_.get())
        throw runtime_error("[PopulationConfigGenerator_IslandBase] Population size trajectory not set.");

    // initialize migration_rate_trajectories_

    TrajectoryPtr zero = TrajectoryPtr(new Trajectory_Constant("id_dummy_zero", 0.));
    migration_rate_trajectories_.resize(population_count_);

    for (MigrationRateTrajectories::iterator it=migration_rate_trajectories_.begin(); it!=migration_rate_trajectories_.end(); ++it)
    {
        it->resize(population_count_);
        for (TrajectoryPtrs::iterator jt=it->begin(); jt!=it->end(); ++jt)
            *jt = zero;
    }

    if (migration_rate_trajectory_default_.get())
        initialize_default_trajectories(); // initialization in derived class

    for (MigrationRateTrajectoryInfos::const_iterator info=migration_rate_trajectory_infos_.begin();
        info!=migration_rate_trajectory_infos_.end(); ++info)
    {
        if (info->population_index_from >= population_count_ ||
            info->population_index_to >= population_count_)
            throw runtime_error(("[PopulationConfigGenerator_IslandBase] Invalid population index:\n    "
                                + info->configuration()).c_str());

        migration_rate_trajectories_[info->population_index_to][info->population_index_from] = info->trajectory;
    }
    
    // check pointers

    for (MigrationRateTrajectories::iterator it=migration_rate_trajectories_.begin(); it!=migration_rate_trajectories_.end(); ++it)
    for (TrajectoryPtrs::iterator jt=it->begin(); jt!=it->end(); ++jt)
    if (!jt->get())
    {
        ostringstream oss;
        oss << "[PopulationConfigGenerator_LinearSteppingStone] Null trajectory from population " 
            << jt - it->begin() << " to population " << it - migration_rate_trajectories_.begin() << endl;
        throw runtime_error(oss.str().c_str());
    }
}


//
// PopulationConfigGenerator_LinearSteppingStone
//


PopulationConfigGenerator_LinearSteppingStone::PopulationConfigGenerator_LinearSteppingStone(const string& id)
:   PopulationConfigGenerator_IslandBase(id)
{}


void PopulationConfigGenerator_LinearSteppingStone::initialize_default_trajectories()
{
    // default migration between adjacent populations only

    for (size_t to=0; to<population_count_; ++to)
        for (size_t from=0; from<population_count_; ++from)
            if (to == from+1 || to+1 == from)
                migration_rate_trajectories_[to][from] = migration_rate_trajectory_default_;
}


//
// PopulationConfigGenerator_Island
//


PopulationConfigGenerator_Island::PopulationConfigGenerator_Island(const string& id)
:   PopulationConfigGenerator_IslandBase(id)
{}


void PopulationConfigGenerator_Island::initialize_default_trajectories()
{
    // default migration between any two populations

    for (size_t to=0; to<population_count_; ++to)
        for (size_t from=0; from<population_count_; ++from)
            if (to != from)
                migration_rate_trajectories_[to][from] = migration_rate_trajectory_default_;
}


  
//
// Note for future PCG implementations:
//   population_configs() has the requirement that population_configs(0) always returns the same value.
//   This is to allow other modules to obtain and rely upon information about the initial populations
//   via initialize(simconfig), e.g. VI_SingleLocusHardyWeinberg.  As a consequence,
//   PCG implementations using random population sizes must cache the result from population_configs(0)
//   and return the same result on future calls.  (This is not an issue for PCG implementations with
//   deterministic population sizes.)
//
    

