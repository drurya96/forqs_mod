//
// subsample_population.cpp
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


#include "Population_Organisms.hpp"
#include "Random.hpp"
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
namespace bfs = boost::filesystem;


struct Config
{
    bfs::path filename;
    size_t sampleSize;
    size_t replicateCount;
    bfs::path outputDirectory;

    Config()
    :   sampleSize(0), replicateCount(0)
    {}
};


void subsamplePopulation(const Config& config)
{
    Population_Organisms p;
    bfs::ifstream is(config.filename);

    cout << "Reading population data.\n";
    is >> p;
    is.close();
    if (p.organisms().empty())
        throw runtime_error("Error reading population data.");

    for (size_t i=0; i<config.replicateCount; i++)
    {
        shared_ptr<Population> subsample = p.randomSubsample(config.sampleSize);

        ostringstream filename;
        filename << "subsample_" << config.sampleSize << "_" << i << ".txt";
        bfs::ofstream os(config.outputDirectory / filename.str());

        os << *subsample;
        os.close();
    }
}


Config parseCommandLine(int argc, char* argv[])
{
    if (argc != 5)
    {
        cout << "Usage: subsample_population <filename> <sample_size> <replicate_count> <outputdir>\n";
        cout << "\n";
        cout << "Darren Kessner\n";
        cout << "John Novembre Lab, UCLA\n";
        throw runtime_error("");
    }

    Config config;
    config.filename = argv[1];
    config.sampleSize = atoi(argv[2]);
    config.replicateCount = atoi(argv[3]);
    config.outputDirectory = argv[4];

    if (!bfs::exists(config.filename))
    {
        ostringstream oss;
        oss << "File not found: " << config.filename;
        throw runtime_error(oss.str());
    }

    if (bfs::exists(config.outputDirectory))
    {
        ostringstream oss;
        oss << "File/directory already exists: " << config.outputDirectory;
        throw runtime_error(oss.str());
    }

    bfs::create_directories(config.outputDirectory);

    return config;
}


int main(int argc, char* argv[])
{
    try
    {
        Config config = parseCommandLine(argc, argv); 
        subsamplePopulation(config);
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


