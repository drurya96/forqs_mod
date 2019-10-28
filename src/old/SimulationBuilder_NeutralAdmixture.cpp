//
// SimulationBuilder_NeutralAdmixture.cpp
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


#include "SimulationBuilder_NeutralAdmixture.hpp"
#include "Random.hpp"
#include "boost/filesystem.hpp"
#include <iostream>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;


class Reporter_Log : public Reporter // simple reporter for testing
{
    public:

    Reporter_Log(const string& output_directory)
    :   Configurable("reporter_log"), outdir_(output_directory)
    {
        /*
        os_log_.open(outdir_ / "log.txt");
        if (!os_log_)
            throw runtime_error("[Reporter_Log] Unable to open log.");
            */
    }

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDatas& population_datas,
                        bool is_final_generation)
    {
        if (!os_log_)
        {
            os_log_.open(outdir_ / "log.txt");
            if (!os_log_)
                throw runtime_error("[Reporter_Log] Unable to open log.");
        }

        if (is_final_generation)
            os_log_ << generation_number << endl;
        else
        {
            // output the last generation

            for (size_t i=0; i<populations.size(); ++i)
            {
                ostringstream filename;
                filename << "pop" << i << ".txt"; 
                bfs::ofstream os_pop(outdir_ / filename.str());
                if (!os_pop)
                    throw runtime_error(("[Reporter_Log] Unable to open " + filename.str()).c_str());
                os_pop << *populations[i];
            }
        }
    }

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_Log";}

    private:

    bfs::path outdir_;
    bfs::ofstream os_log_;
};


SimulationBuilder_NeutralAdmixture::Config::Config(const Parameters& parameters)
{
    seed = parameters.value<size_t>("seed", 0);
    population_config_filename = parameters.value<string>("popconfig", "");
    genetic_map_list_filename = parameters.value<string>("genetic_map_list", "");
    output_directory = parameters.value<string>("outdir", "");
}


SimulationBuilder_NeutralAdmixture::SimulationBuilder_NeutralAdmixture(const SimulationBuilder_NeutralAdmixture::Config& config)
:   config_(config)
{}


void SimulationBuilder_NeutralAdmixture::usage() const
{
    cout << "Usage:  simrecomb neutral_admixture [parameter_name=value] ...\n";
    cout << "        simrecomb na [parameter_name=value] ...\n";
    cout << "        simrecomb na config=config_filename ...\n";
    cout << endl;
    cout << "Required parameters:\n";
    cout << "  outdir=<output_directory>\n";
    cout << "  popconfig=<population_config_filename>\n";
    cout << "  genetic_map_list=<genetic_map_list_filename>\n";
    cout << endl;
    cout << "Optional parameters:\n";
    cout << "  config=<config_filename>\n";
    cout << "  seed=<value>\n";
    cout << endl;
}


SimulatorConfigPtr SimulationBuilder_NeutralAdmixture::create_simulator_config() const
{
    // create output directory

    /*
    if (config_.output_directory.empty())
        throw runtime_error("[SimulationBuilder_NeutralAdmixture] No output directory specified (outdir=value).");

    if (bfs::exists(config_.output_directory))
        throw runtime_error(("[SimulationBuilder_NeutralAdmixture] Output directory exists: " + config_.output_directory).c_str());

    bfs::create_directories(config_.output_directory);
    */

    // read configuration files

    if (config_.population_config_filename.empty())
        throw runtime_error("[SimulationBuilder_NeutralAdmixture] No population config file specified (popconfig=value).");

    if (!bfs::exists(config_.population_config_filename))
        throw runtime_error(("[SimulationBuilder_NeutralAdmixture] Population config file not found: " + config_.population_config_filename).c_str());

    SimulatorConfigPtr simulator_config(new SimulatorConfig("simconfig"));

    cout << "[SimulationBuilder_NeutralAdmixture] Reading population configuration file " << config_.population_config_filename << endl;
    simulator_config->population_config_generator = PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_File("pcg_file", config_.population_config_filename));

    cout << "[SimulationBuilder_NeutralAdmixture] Reading genetic map list " << config_.genetic_map_list_filename << endl;
    vector<string> genetic_map_filenames;
    bfs::ifstream is_genetic_map_list(config_.genetic_map_list_filename);
    copy(istream_iterator<string>(is_genetic_map_list), istream_iterator<string>(), back_inserter(genetic_map_filenames));
    is_genetic_map_list.close();

    // initialize simulator

    cout << "seed: " << config_.seed << endl;
    cout << "genetic maps:\n";
    copy(genetic_map_filenames.begin(), genetic_map_filenames.end(), ostream_iterator<string>(cout, "\n"));
    cout << endl;

    // initialize recombination maps

    cout << "[SimulationBuilder_NeutralAdmixture] Initializing recombination maps.\n";

    //Random::seed(config_.seed);

    simulator_config->recombination_position_generator = RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_RecombinationMap("dummy_id", 
                                                                genetic_map_filenames));

    cout << "output_directory: " << config_.output_directory << endl;

    simulator_config->reporters.push_back(ReporterPtr(new Reporter_Log(config_.output_directory)));
    simulator_config->reporters.push_back(ReporterPtr(new Reporter_Timer(config_.output_directory)));

    simulator_config->seed = config_.seed;
    simulator_config->output_directory = config_.output_directory;

    return simulator_config;
}

