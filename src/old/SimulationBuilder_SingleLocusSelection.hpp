//
// SimulationBuilder_SingleLocusSelection.hpp
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


#ifndef _SIMULATIONBUILDER_SINGLELOCUSSELECTION_HPP_
#define _SIMULATIONBUILDER_SINGLELOCUSSELECTION_HPP_


#include "SimulationBuilder.hpp"
#include "Parameters.hpp"


class SimulationBuilder_SingleLocusSelection : public SimulationBuilder
{
    public:

    struct Config
    {
        int seed;                               // "seed"
        std::string output_directory;           // "outdir"
        size_t population_count;                // "popcount"
        size_t population_size;                 // "popsize"
        size_t generation_count;                // "gencount"
        double initial_allele_frequency;        // "allelefreq"
        std::vector<double> w;                  // "w0", "w1", "w2" (relative fitnesses for genotype in {0,1,2})
        bool verbose;                           // "verbose" (for debugging)

        Config(const Parameters& parameters = Parameters()); // allows auto conversion: Parameters->Config
    };

    SimulationBuilder_SingleLocusSelection(const Config& config = Config());

    virtual void usage() const;
    virtual SimulatorConfigPtr create_simulator_config() const;

    private:

    Config config_;
};


std::ostream& operator<<(std::ostream& os, const SimulationBuilder_SingleLocusSelection::Config& config);


#endif //  _SIMULATIONBUILDER_SINGLELOCUSSELECTION_HPP_

