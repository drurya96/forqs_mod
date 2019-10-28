//
// PopulationConfigGeneratorExperimental.hpp
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


#ifndef _POPULATIONCONFIGGENERATOREXPERIMENTAL_HPP_
#define _POPULATIONCONFIGGENERATOREXPERIMENTAL_HPP_


#include "PopulationConfigGenerator.hpp"


//
// PopulationConfigGenerator_TurnerExperiment
//

class PopulationConfigGenerator_TurnerExperiment : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_TurnerExperiment(const std::string& id);

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "PopulationConfigGenerator_TurnerExperiment";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);

    private:

    size_t population_size_founders_;
    size_t population_size_neutral_;
    size_t population_size_replicate_;
    size_t generation_count_neutral_;
    size_t generation_count_selected_;
    std::string qtid_;
    double proportion_selected_;
    bool report_hidden_;
};


#endif // _POPULATIONCONFIGGENERATOREXPERIMENTAL_HPP_

