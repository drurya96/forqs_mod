//
// MutationGenerator.hpp
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


#ifndef _MUTATIONGENERATOR_HPP_
#define _MUTATIONGENERATOR_HPP_


#include "Locus.hpp"
#include "Population.hpp"
#include "shared_ptr.hpp"


//
// MutationGenerator
//

///
/// interface class
///

class MutationGenerator : public Configurable
{
    public:

    struct MutationInfo
    {
        size_t individual_index;
        Locus locus;
        int which;              // 0 == chromosome_pair.first, 1 == chromosome_pair.second
        unsigned int value;     // new value to be returned by VariantIndicator

        MutationInfo()
        :   individual_index(0), locus("id_dummy"), which(0), value(0)
        {}
    };

    typedef std::vector<MutationInfo> MutationInfos;

    MutationGenerator(const std::string& id) : Configurable(id) {}

    virtual MutationInfos generate_mutations(const Population& population,
                                             size_t generation_index,
                                             size_t population_index) const {return MutationInfos();}
    // Configurable interface

    virtual std::string class_name() const {return "MutationGenerator";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}

    virtual ~MutationGenerator() {} 
};


typedef boost::shared_ptr<MutationGenerator> MutationGeneratorPtr;


std::ostream& operator<<(std::ostream& os, const MutationGenerator::MutationInfo& mutation_info);


#endif // _MUTATIONGENERATOR_HPP_

