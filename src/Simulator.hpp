//
// Simulator.hpp
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


#ifndef _SIMULATOR_HPP_
#define _SIMULATOR_HPP_


#include "PopulationConfigGenerator.hpp"
#include "PopulationData.hpp"
#include "MutationGenerator.hpp"
#include "QuantitativeTrait.hpp"
#include "Reporter.hpp"
#include "VariantIndicator.hpp"
#include "FitnessFunctionImplementation.hpp"
#include <vector>
#include <string>
#include <iostream>


// forward declarations to make SimulatorConfig accessible to VariantIndicator and Reporter
class VariantIndicator;
typedef boost::shared_ptr<VariantIndicator> VariantIndicatorPtr;
class Reporter;
typedef boost::shared_ptr<Reporter> ReporterPtr;
typedef std::vector<ReporterPtr> ReporterPtrs;


///
/// \defgroup SimulatorConfig SimulatorConfig
///
/// main simulator configuration
///


///
/// main data structure holding global simulation parameters and references to top-level modules
///
/// Global simulation parameters:
/// parameter | default | notes
/// ----------|---------|-------------
/// output_directory = \<string\> | none | required
/// seed = \<float\> | 0 | optional
/// write_popconfig = \<int\> | 0 (= don't write) | optional
///
/// References to top-level modules:
/// parameter | default | notes
/// ----------|---------|-------------
/// population_config_generator = \<id\> | none | required
/// recombination_position_generator = \<id\> | RecombinationPositionGenerator_Trivial | optional
/// variant_indicator = \<id\> | VariantIndicator_Trivial | optional
/// quantitative_trait = \<id\> | none | optional, multiple ok
/// mutation_generator = \<id\> | none | optional
/// reporter = \<id\> | none | optional, multiple ok
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup SimulatorConfig
///

struct SimulatorConfig : public Configurable
{
    unsigned int seed;
    std::string output_directory;
    bool write_popconfig;
    bool write_vi;
    bool use_random_seed;
	bool recombination_tracker = false;

    PopulationConfigGeneratorPtr population_config_generator;
    RecombinationPositionGeneratorPtrs recombination_position_generators;
	RecombinationPositionGeneratorPtrsArray recombination_position_generators_array;
    VariantIndicatorPtr variant_indicator;
    QuantitativeTraitPtrs quantitative_traits;
    MutationGeneratorPtr mutation_generator;
    ReporterPtrs reporters;
	Configurable::Registry registry;

    SimulatorConfig(const std::string& id = "dummy");

    // Configurable interface

    virtual std::string class_name() const {return "SimulatorConfig";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;
};

typedef boost::shared_ptr<SimulatorConfig> SimulatorConfigPtr;

class Simulator
{
    public:

    Simulator(SimulatorConfig& config, const Parameters& parameters);
    void simulate_single_generation();
    void simulate_all();
    void update_final();

    private:

    SimulatorConfig config_;
    Genotyper genotyper_;

    size_t current_generation_index_;
    PopulationPtrsPtr current_populations_;
    PopulationDataPtrsPtr current_population_datas_;
    size_t update_step_;

    bfs::ofstream os_popconfigs_;
};


typedef boost::shared_ptr<Simulator> SimulatorPtr;

void set_recombination_generators(FitnessFunction_Recombination& ff, RecombinationPositionGeneratorPtrsArray& rpg_vector);

#endif //  _SIMULATOR_HPP_

