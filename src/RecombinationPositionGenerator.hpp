//
// RecombinationPositionGenerator.hpp
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


#ifndef _RECOMBINATIONPOSITIONGENERATOR_HPP_
#define _RECOMBINATIONPOSITIONGENERATOR_HPP_


#include "Configurable.hpp"
#include "shared_ptr.hpp"
#include <vector>


//
// RecombinationPositionGenerator 
//

///
/// interface class
///

class RecombinationPositionGenerator : public Configurable
{
    public:

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index = 0) const = 0;
    virtual ~RecombinationPositionGenerator(){}

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    protected:

    RecombinationPositionGenerator(const std::string& id) : Configurable(id) {}
};


typedef boost::shared_ptr<RecombinationPositionGenerator> RecombinationPositionGeneratorPtr;
typedef std::vector<RecombinationPositionGeneratorPtr> RecombinationPositionGeneratorPtrs;

typedef RecombinationPositionGeneratorPtrsArray* RecombinationPositionGeneratorPtrs; // A.D. 10/29/2019


#endif // _RECOMBINATIONPOSITIONGENERATOR_HPP_

