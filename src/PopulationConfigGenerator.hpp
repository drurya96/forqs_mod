//
// PopulationConfigGenerator.hpp
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


#ifndef _POPULATIONCONFIGGENERATOR_HPP_
#define _POPULATIONCONFIGGENERATOR_HPP_


#include "Configurable.hpp"
#include "Population.hpp"
#include "PopulationData.hpp"
#include "shared_ptr.hpp"


//
// PopulationConfigGenerator 
//

///
/// interface class
///


class PopulationConfigGenerator : public Configurable
{
    public:

    size_t generation_count() const {return generation_count_;}
    size_t population_count() const {return population_count_;}
    unsigned int id_offset_step() const {return id_offset_step_;}
    size_t chromosome_pair_count() const {return chromosome_pair_count_;}
    const std::vector<unsigned int>& chromosome_lengths() const {return chromosome_lengths_;}
    const std::vector<std::string>& fitness_functions() const {return fitness_functions_;}

    virtual Population::Configs population_configs(size_t generation_index, 
                                                   const PopulationDataPtrs& population_datas) const = 0;

    unsigned int min_unused_id() const;

    // Configurable interface: base implementation

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    virtual ~PopulationConfigGenerator() {}

    protected:

    PopulationConfigGenerator(const std::string& id);

    size_t generation_count_;
    size_t population_count_;
    unsigned int id_offset_step_;
    size_t chromosome_pair_count_;
    std::vector<unsigned int> chromosome_lengths_;
    std::vector<std::string> fitness_functions_;
};


typedef shared_ptr<PopulationConfigGenerator> PopulationConfigGeneratorPtr;


#endif // _POPULATIONCONFIGGENERATOR_HPP_

