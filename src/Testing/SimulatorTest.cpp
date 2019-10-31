//
// SimulatorTest.cpp
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


#include "Simulator.hpp"
#include "PopulationConfigGeneratorImplementation.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "VariantIndicatorImplementation.hpp"
#include "QuantitativeTraitImplementation.hpp"
#include "FitnessFunctionImplementation.hpp"
#include "ReporterImplementation.hpp"
#include "unit.hpp"
#include <iostream>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_Configurable_SimulatorConfig()
{
    if (os_) *os_ << "test_Configurable_SimulatorConfig()\n";

    Configurable::Registry registry;
    registry["my_pcg"] = PopulationConfigGeneratorPtr(new PopulationConfigGenerator_ConstantSize("my_pcg"));
    registry["my_rpg"] = RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Trivial("my_rpg"));
    registry["my_rpg_2"] = RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Trivial("my_rpg_2"));
    registry["my_variant_indicator"] = VariantIndicatorPtr(new VariantIndicator_Trivial("my_variant_indicator"));
    registry["my_qt_1"] = QuantitativeTraitPtr(new QuantitativeTrait_SingleLocusFitness("my_qt_1"));
    registry["my_qt_2"] = QuantitativeTraitPtr(new QuantitativeTrait_SingleLocusFitness("my_qt_2"));
    registry["my_ff"] = QuantitativeTraitPtr(new FitnessFunction_Trivial("my_ff"));
    registry["my_mutation_generator"] = MutationGeneratorPtr(new MutationGenerator("my_mutation_generator"));
    registry["my_reporter_1"] = ReporterPtr(new Reporter_AlleleFrequencies("my_reporter_1"));

    Parameters parameters_in;
    parameters_in.insert_name_value("population_config_generator", "my_pcg");
    parameters_in.insert_name_value("recombination_position_generator", "my_rpg");
    parameters_in.insert_name_value("recombination_position_generator", "my_rpg_2");
    parameters_in.insert_name_value("variant_indicator", "my_variant_indicator");
    parameters_in.insert_name_value("quantitative_trait", "my_qt_1");
    parameters_in.insert_name_value("quantitative_trait", "my_qt_2");
    parameters_in.insert_name_value("mutation_generator", "my_mutation_generator");
    parameters_in.insert_name_value("reporter", "my_reporter_1");
    parameters_in.insert_name_value("seed", 1234);
    parameters_in.insert_name_value("output_directory", "blah");
    parameters_in.insert_name_value("write_popconfig", true);
    parameters_in.insert_name_value("write_vi", true);

    SimulatorConfig config("dummy_id");
    config.configure(parameters_in, registry);

    Parameters parameters_out = config.parameters();

    if (os_)
    {
        *os_ << "parameters_in:\n" << parameters_in << endl
             << "parameters_out:\n" << parameters_out << endl
             << "configuration:\n\n";
        set<string> ids_written;
        config.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters_in == parameters_out);
}


void test()
{
    test_Configurable_SimulatorConfig();
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


