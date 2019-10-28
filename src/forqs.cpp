//
// forqs.cpp
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


#include "SimulationBuilder_Generic.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;


void parse_command_line(int argc, char* argv[], string& filename, Parameters& parameters)
{
    ostringstream usage;

    usage << "Forward-in-time simulation of Recombination, Quantitative traits, and Selection\n"
          << "Created by Darren Kessner with John Novembre at UCLA\n"
          << "Copyright (c) 2013 Regents of the University of California\n"
          << endl
          << "Usage: forqs config_file [parameters]\n";

    if (argc < 2)
        throw runtime_error(usage.str().c_str());

    filename = argv[1];

    for (int i=2; i<argc; ++i)
        parameters.parse(argv[i]); // parse any other parameters
}


int main(int argc, char* argv[])
{
    cout << "forqs\n\n";

    try
    {
        string filename;
        Parameters parameters;
        parse_command_line(argc, argv, filename, parameters);

        SimulationBuilder_Generic builder(filename, parameters);
        SimulatorConfigPtr simconfig = builder.create_simulator_config();

        Simulator simulator(*simconfig);
        simulator.simulate_all();
        simulator.update_final();

        builder.write_new_seed();

        return 0;
    }
    catch(exception& e)
    {
        cerr << endl << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << endl << "Caught unknown exception.\n";
        return 1;
    }
}


