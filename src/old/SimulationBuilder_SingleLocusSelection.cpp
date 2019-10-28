//
// SimulationBuilder_SingleLocusSelection.cpp
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


#include "SimulationBuilder_SingleLocusSelection.hpp"
#include "Random.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/lambda/lambda.hpp"
#include <iostream>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;
using namespace boost::lambda;


SimulationBuilder_SingleLocusSelection::Config::Config(const Parameters& parameters)
:   w(3, 1.0)
{
    seed = parameters.value<int>("seed", 1);
    output_directory = parameters.value<string>("outdir", "");

    population_count = parameters.value<size_t>("popcount", 1);
    population_size = parameters.value<size_t>("popsize", 0);
    generation_count = parameters.value<size_t>("gencount", 0);
    initial_allele_frequency = parameters.value<double>("allelefreq", 0);

    w[0] = parameters.value<double>("w0", 1.0);
    w[1] = parameters.value<double>("w1", 1.0);
    w[2] = parameters.value<double>("w2", 1.0);

    verbose = parameters.value<bool>("verbose", false);
}


std::ostream& operator<<(std::ostream& os, const SimulationBuilder_SingleLocusSelection::Config& config)
{
    os << "seed = " << config.seed << endl;
    os << "outdir = " << config.output_directory << endl;
    os << "popcount = " << config.population_count << endl;
    os << "popsize = " << config.population_size << endl;
    os << "gencount = " << config.generation_count << endl;
    os << "allelefreq = " << config.initial_allele_frequency << endl;
    os << "w0 = " << config.w[0] << endl;
    os << "w1 = " << config.w[1] << endl;
    os << "w2 = " << config.w[2] << endl;
    if (config.verbose) os << "verbose" << endl;
    return os;
}


SimulationBuilder_SingleLocusSelection::SimulationBuilder_SingleLocusSelection(const Config& config)
:   config_(config)
{}


void SimulationBuilder_SingleLocusSelection::usage() const
{   
    cout << "Usage:  simrecomb single_locus_selection [parameter_name=value] ...\n";
    cout << "        simrecomb sls [parameter_name=value] ...\n";
    cout << "        simrecomb sls config=config_filename ...\n";
    cout << endl;
    cout << "Required parameters:\n";
    cout << "  outdir=<output_directory>\n";
    cout << "  popsize=<population_size>\n";
    cout << "  gencount=<generation_count>\n";
    cout << endl;
    cout << "Optional parameters:\n";
    cout << "  config=<config_filename>\n";
    cout << "  seed=<value>\n";
    cout << "  popcount=<population_count>              (default: 1)\n";
    cout << "  allelefreq=<initial_allele_frequency>    (default: allelefreq=0)\n";
    cout << "  w0=<relative_fitness_genotype_0>         (default: w0=1)\n";
    cout << "  w1=<relative_fitness_genotype_1>         (default: w1=1)\n";
    cout << "  w2=<relative_fitness_genotype_2>         (default: w2=1)\n";
    cout << "  verbose\n"; 
    cout << endl;
}


SimulatorConfigPtr SimulationBuilder_SingleLocusSelection::create_simulator_config() const
{
    // check parameters

    if (config_.output_directory.empty())
        throw runtime_error("[SimulationBuilder_SingleLocusSelection] No output directory specified (outdir=value).");

    if (bfs::exists(config_.output_directory))
        throw runtime_error(("[SimulationBuilder_SingleLocusSelection] Output directory exists: " + config_.output_directory).c_str());

    if (config_.population_size == 0)
        throw runtime_error("[SimulationBuilder_SingleLocusSelection] Population size not specified (popsize=value).");

    if (config_.generation_count == 0)
        throw runtime_error("[SimulationBuilder_SingleLocusSelection] Generation count not specified (gencount=value).");

    //bfs::create_directories(config_.output_directory);

    bfs::ofstream os_config(bfs::path(config_.output_directory) / "config.txt");
    os_config << config_ << endl;
    os_config.close();

    cout << config_ << endl;

    // initialize

    SimulatorConfigPtr simulator_config(new SimulatorConfig("simconfig"));

    /*
    PopulationConfigGenerator_ConstantSize::Config constant_size_config;
    constant_size_config.chromosome_pair_count = 1;
    constant_size_config.population_size = config_.population_size;
    constant_size_config.population_count = config_.population_count;
    constant_size_config.generation_count = config_.generation_count;
    */

    simulator_config->population_config_generator = PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_ConstantSize("popconfig_generator"));

    Parameters parameters_pcg;
    parameters_pcg.insert_name_value("chromosome_pair_count", 1);
    parameters_pcg.insert_name_value("population_size", config_.population_size);
    parameters_pcg.insert_name_value("population_count", config_.population_count);
    parameters_pcg.insert_name_value("generation_count", config_.generation_count);
    parameters_pcg.insert_name_value("chromosome_lengths", "1000000");
    
    Configurable::Registry registry;
    simulator_config->population_config_generator->configure(parameters_pcg, registry);

    Locus locus("locus", 0, 500000); // hardcoded locus
    string qtid = "my_qtid"; // hardcoded QT id

    simulator_config->recombination_position_generator = RecombinationPositionGeneratorPtr(
        new RecombinationPositionGenerator_SingleCrossover("rpg_uniform"));

    simulator_config->variant_indicator = VariantIndicatorPtr(new VariantIndicator_SingleLocusHardyWeinberg(
        "my_variant_indicator", locus, config_.initial_allele_frequency));
    
    QuantitativeTraitPtr qt(new QuantitativeTrait_SingleLocusFitness(qtid, locus, config_.w));
    simulator_config->quantitative_traits.push_back(qt);

    FitnessFunctionPtr ff(new FitnessFunction_Identity("ff_identity", qtid));
    simulator_config->fitness_function = ff;

    if (config_.verbose)
        simulator_config->reporters.push_back(ReporterPtr(new Reporter_Population("reporter_population")));
    //simulator_config->reporters.push_back(ReporterPtr(new Reporter_AlleleFrequencies("reporter_allele_frequencies", locus)));
    simulator_config->reporters.push_back(ReporterPtr(new Reporter_MeanFitnesses("reporter_mean_fitnesses")));
    simulator_config->reporters.push_back(ReporterPtr(new Reporter_Timer("reporter_timer")));

    simulator_config->reporters.push_back(ReporterPtr(
        new Reporter_DeterministicTrajectories("reporter_deterministic_trajectories",
                                               config_.initial_allele_frequency,
                                               config_.w)));

    simulator_config->seed = config_.seed;
    simulator_config->output_directory = config_.output_directory;

    return simulator_config;    
}


