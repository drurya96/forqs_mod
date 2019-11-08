//
// PopulationConfigGeneratorImplementation.hpp
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


#ifndef _POPULATIONCONFIGGENERATORIMPLEMENTATION_HPP_
#define _POPULATIONCONFIGGENERATORIMPLEMENTATION_HPP_


#include "PopulationConfigGenerator.hpp"
#include "Trajectory.hpp"

#include "Configurable.hpp"

///
/// \defgroup PopulationConfigGenerators PopulationConfigGenerators
///
/// classes for specifying demography
///
/// PopulationConfigGenerator implementations generate a population configuration
/// (Population::Config) for each generation.  The population configuration includes
///   - population sizes
///   - mating distributions (joint distribution of parent populations from which to initialize offspring populations)
///


//
// PopulationConfigGenerator_File
//

///
/// reads full specification of each Population::Config from a file
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// filename = \<population_config_filename\> | none | required
///
/// Example: [example_neutral_admixture.txt](../../examples/example_neutral_admixture.txt)
///
/// \ingroup PopulationConfigGenerators
///

class PopulationConfigGenerator_File : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_File(const std::string& id, const std::string& filename = "");
    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "PopulationConfigGenerator_File";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::string filename_;
    std::vector<Population::Configs> population_configs_; 

    void read_file();
};


//
// PopulationConfigGenerator_ConstantSize
//

///
/// generates isolated populations of constant size
/// 
/// PopulationConfigGenerator general parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// generation_count = \<int\> | none | required
/// population_count = \<int\> | 1 | optional
/// id_offset_step = \<int\> | 0 | optional
/// chromosome_pair_count = \<int\> | 1 | optional
/// chromosome_lengths = \<int\> [\<int\> ... ] | none | optional (must match chromosome_pair_count)
///
/// PopulationConfigGenerator_ConstantSize specific parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// population_size = \<int\> | none | required
///
/// Note: the populations are isolated (i.e. no migration)
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup PopulationConfigGenerators
///

class PopulationConfigGenerator_ConstantSize : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_ConstantSize(const std::string& id)
    :   PopulationConfigGenerator(id)
    {}

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "PopulationConfigGenerator_ConstantSize";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:
    size_t population_size_;
};


//
// MigrationRateTrajectoryInfo
//


struct MigrationRateTrajectoryInfo
{
    size_t population_index_from;
    size_t population_index_to;
    TrajectoryPtr trajectory;
    
    MigrationRateTrajectoryInfo() : population_index_from(0), population_index_to(0) {}
    MigrationRateTrajectoryInfo(const std::string& configuration, const Configurable::Registry& registry);
	std::string configuration() const;
};


typedef std::vector<MigrationRateTrajectoryInfo> MigrationRateTrajectoryInfos;


//
// PopulationConfigGenerator_IslandBase 
//
// This base implmentation class provides common functionality for
// PCG_LinearSteppingStone and PCG_Island, including parametrization, internal
// data structures, and population_configs() implementation.
//
// Note that the only difference between these models is the interpretation of
// the migration_rate_default parameter; specialization is handled by the
// virtual function initialize_default_trajectories().
//

class PopulationConfigGenerator_IslandBase : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_IslandBase(const std::string& id);

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    protected:

    TrajectoryPtr population_size_trajectory_;
    TrajectoryPtr migration_rate_trajectory_default_;
    MigrationRateTrajectoryInfos migration_rate_trajectory_infos_;

    typedef std::vector<TrajectoryPtrs> MigrationRateTrajectories;
    MigrationRateTrajectories migration_rate_trajectories_;

    virtual void initialize_default_trajectories() = 0;

    private:

    void initialize_trajectories();
};


//
// PopulationConfigGenerator_LinearSteppingStone
//

///
/// generates population configurations for a linear stepping stone model,
/// with specified trajectories for population sizes and migration rates
/// 
/// PopulationConfigGenerator general parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// generation_count = \<int\> | none | required
/// population_count = \<int\> | 1 | optional 
/// id_offset_step = \<int\> | 0 | optional
/// chromosome_pair_count = \<int\> | 1 | optional
/// chromosome_lengths = \<int\> [\<int\> ... ] | none | optional (must match chromosome_pair_count)
///
/// PopulationConfigGenerator_LinearSteppingStone specific parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// population_size = \<id_trajectory\> | none | required
/// migration_rate_default = \<id_trajectory\> | none | optional
/// migration_rate:from:to = \<id_trajectory\> \<int_population_from\> \<int_population_to\> | none | optional (1 == first population)
///
/// Notes: 
/// - if specified, migration_rate_default is the default migration rate between adjacent populations
/// - id_offset_step is used only for assigning ids to initial populations (useful only in special cases)
///
/// Example: [example_stepping_stone.txt](../../examples/example_stepping_stone.txt)
///
/// \ingroup PopulationConfigGenerators
///


class PopulationConfigGenerator_LinearSteppingStone : public PopulationConfigGenerator_IslandBase
{
    public:

    PopulationConfigGenerator_LinearSteppingStone(const std::string& id);
    virtual std::string class_name() const {return "PopulationConfigGenerator_LinearSteppingStone";}

    protected:

    virtual void initialize_default_trajectories();
};


//
// PopulationConfigGenerator_Island
//

///
/// generates population configurations for an island model,
/// with specified trajectories for population sizes and migration rates
///
/// PopulationConfigGenerator general parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// generation_count = \<int\> | none | required
/// population_count = \<int\> | 1 | optional
/// id_offset_step = \<int\> | 0 | optional
/// chromosome_pair_count = \<int\> | 1 | optional
/// chromosome_lengths = \<int\> [\<int\> ... ] | none | optional (must match chromosome_pair_count)
///
/// PopulationConfigGenerator_Island specific parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// population_size = \<id_trajectory\> | none | required
/// migration_rate_default = \<id_trajectory\> | none | optional
/// migration_rate:from:to = \<id_trajectory\> \<int_population_from\> \<int_population_to\> | none | optional (1 == first population)
///
/// Notes: 
/// - if specified, migration_rate_default is the default migration rate between any two populations
/// - id_offset_step is used only for assigning ids to initial populations (useful only in special cases)
///
/// Example: [example_stepping_stone_using_island.txt](../../examples/example_stepping_stone_using_island.txt)
///
/// \ingroup PopulationConfigGenerators
///

class PopulationConfigGenerator_Island : public PopulationConfigGenerator_IslandBase
{
    public:

    PopulationConfigGenerator_Island(const std::string& id);
    virtual std::string class_name() const {return "PopulationConfigGenerator_Island";}

    protected:

    virtual void initialize_default_trajectories();
};


#endif // _POPULATIONCONFIGGENERATORIMPLEMENTATION_HPP_

