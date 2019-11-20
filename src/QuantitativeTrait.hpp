//
// QuantitativeTrait.hpp
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


#ifndef _QUANTITATIVETRAIT_HPP_
#define _QUANTITATIVETRAIT_HPP_


#include "Locus.hpp"
#include "Configurable.hpp"
#include "PopulationData.hpp"
#include "shared_ptr.hpp"

typedef boost::shared_ptr<Population> PopulationPtr;
typedef std::vector<PopulationPtr> PopulationPtrs;

//
// QuantitativeTrait
//

///
/// interface class
///

class QuantitativeTrait : public Configurable
{
    public:

    // return set of loci (the QTLs contributing to the trait)
    const Loci& loci() const {return loci_;}
    // calculate trait values for a single population using genotypes
    virtual void calculate_trait_values(const PopulationData& population_data, const Population& population) const;

    // calculate trait values for all populations
    //
    // note: default implementation (simple iteration over populations)
    //       will work for most traits; however, some traits may need information
    //       about previously calculated trait values in all populations
    //       (e.g. threshold for truncation fitness)
    virtual void calculate_trait_values(const PopulationDataPtrs& population_datas, const PopulationPtrs& populations) const;

    virtual ~QuantitativeTrait() {}

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    protected:

    QuantitativeTrait(const std::string& id) : Configurable(id) {}

    Loci loci_;
};


typedef boost::shared_ptr<QuantitativeTrait> QuantitativeTraitPtr;
typedef std::vector<QuantitativeTraitPtr> QuantitativeTraitPtrs;


#endif //  _QUANTITATIVETRAIT_HPP_

