//
// Reporter.cpp
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


#include "Reporter.hpp"
#include "Simulator.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <stdexcept>
#include <numeric>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;


//
// Reporter
//

std::string Reporter::class_name() const
{
    cerr << "[Reporter] Warning: virtual class_name() has not been defined in derived class.\n";
    return "Reporter";
}


Parameters Reporter::parameters() const
{
    cerr << "[Reporter] Warning: virtual parameters() has not been defined in derived class.\n";
    return Parameters();
}


void Reporter::configure(const Parameters& parameters, const Registry& registry)
{
    cerr << "[Reporter] Warning: virtual configure() has not been defined in derived class.\n";
}


void Reporter::initialize(const SimulatorConfig& config)
{
    output_directory_ = config.output_directory;
}


