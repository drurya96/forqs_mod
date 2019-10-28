//
// PopulationConfigGeneratorExperimentalTest.cpp
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


#include "PopulationConfigGeneratorExperimental.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_PopulationConfigGenerator_TurnerExperiment()
{
    if (os_) *os_ << "test_PopulationConfigGenerator_TurnerExperiment()\n";

    PopulationConfigGenerator_TurnerExperiment pcg("pcg");

    Parameters parameters;
    parameters.insert_name_value("population_size_founders", 810);
    parameters.insert_name_value("population_size_neutral", 1000);
    parameters.insert_name_value("population_size_replicate", 331);
    parameters.insert_name_value("generation_count_neutral", 5);
    parameters.insert_name_value("generation_count_selected", 13);
    parameters.insert_name_value("quantitative_trait", "ipi");
    parameters.insert_name_value("proportion_selected", .2);
    parameters.insert_name_value("report_hidden", 1);

    Configurable::Registry registry;
    pcg.configure(parameters, registry);

    Parameters parameters_out = pcg.parameters();
    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test()
{
    test_PopulationConfigGenerator_TurnerExperiment();
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


