//
// MutationGeneratorImplementation.hpp
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


#ifndef _MUTATIONGENERATORIMPLEMENTATION_HPP_
#define _MUTATIONGENERATORIMPLEMENTATION_HPP_


#include "MutationGenerator.hpp"
#include "Trajectory.hpp"


///
/// \defgroup MutationGenerators MutationGenerators
///
/// classes for generating mutations in a population
///


//
// MutationGenerator_SingleLocus
//


///
/// generates mutations at a single locus, with a specified mutation rate
/// (per chromosome, per generation).
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus = \<id\> | none | required
/// mu = \<float\> | none | required
///
/// Note: implementation uses Poisson approximation to obtain the number of
/// mutations in a given generation.
///
/// Example: [example_mutation_selection_additive.txt](../../examples/example_mutation_selection_additive.txt)
///
/// \ingroup MutationGenerators
///


class MutationGenerator_SingleLocus : public MutationGenerator
{
    public:

    MutationGenerator_SingleLocus(const std::string& id,
                                  const Locus& locus = Locus("id_dummy"),
                                  double mu = 0.)
    :   MutationGenerator(id), locus_(locus), mu_(mu)
    {}

    virtual MutationInfos generate_mutations(const Population& population,
                                             size_t generation_index,
                                             size_t population_index) const;
    // Configurable interface

    virtual std::string class_name() const {return "MutationGenerator_SingleLocus";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    Locus locus_;
    double mu_;
};


//
// MutationGenerator_Regions
//


///
/// generates mutations in multiple regions, with specified mutation rate trajectories
/// (per chromosome, per generation).  
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus:length:rate = \<id_locus_start\> \<int_region_length\> \<id_mutation_rate_trajectory\> | none | multiple allowed
/// 
/// Note: implementation uses Poisson approximation to obtain the number of
/// mutations in a given generation.
///
/// Example: [example_neutral_mutation_region.txt](../../examples/example_neutral_mutation_region.txt)
///
/// \ingroup MutationGenerators
///


class MutationGenerator_Regions : public MutationGenerator
{
    public:

    struct RegionInfo
    {
        Locus locus;
        size_t length;
        TrajectoryPtr mutation_rate;

        RegionInfo(const Locus& _locus = Locus("id_dummy"),
                   size_t _length = 0, 
                   TrajectoryPtr _mutation_rate = TrajectoryPtr())
        :   locus(_locus), length(_length), mutation_rate(_mutation_rate)
        {}

        std::string configuration() const;
    };

    typedef std::vector<RegionInfo> RegionInfos;

    MutationGenerator_Regions(const std::string& id,
                              const RegionInfos& region_infos = RegionInfos())
    :   MutationGenerator(id), region_infos_(region_infos)
    {}

    virtual MutationInfos generate_mutations(const Population& population,
                                             size_t generation_index,
                                             size_t population_index) const;
    // Configurable interface

    virtual std::string class_name() const {return "MutationGenerator_Regions";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    RegionInfos region_infos_;
};


#endif // _MUTATIONGENERATORIMPLEMENTATION_HPP_

