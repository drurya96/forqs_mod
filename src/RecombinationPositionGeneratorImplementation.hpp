//
// RecombinationPositionGeneratorImplementation.hpp
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


#ifndef _RECOMBINATIONPOSITIONGENERATORIMPLEMENTATION_HPP_
#define _RECOMBINATIONPOSITIONGENERATORIMPLEMENTATION_HPP_


#include "RecombinationPositionGenerator.hpp"
#include "RecombinationMap.hpp"
#include "Random.hpp"


///
/// \defgroup RecombinationPositionGenerators RecombinationPositionGenerators
///
/// classes generating a list of recombination positions
///
/// RecombinationPositionGenerators are used during the creation of offspring from parents
/// to create recombined chromosomes from the parental chromosomes.
///
/// By convention, zero in the first recombination position indicates that the 
/// recombined chromosome begins with the 2nd chromosome of the parental pair.
/// For example:
///   - empty == full 1st chromosome 
///   - \<0\> == full 2nd chromosome
/// 


//
// RecombinationPositionGenerator_Trivial
//

///
/// trivial implementation: returns one of the two whole chromosomes at random
/// 
/// Parameters: none
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup RecombinationPositionGenerators
///

class RecombinationPositionGenerator_Trivial : public RecombinationPositionGenerator
{
    public:

    RecombinationPositionGenerator_Trivial(const std::string& id) 
    :   RecombinationPositionGenerator(id)
    {}

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "RecombinationPositionGenerator_Trivial";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
};


//
// RecombinationPositionGenerator_SingleCrossover
//

///
/// generates exactly one recombination event, at a position chosen uniformly at random along the chromosome
/// 
/// Parameters: none
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup RecombinationPositionGenerators
///


class RecombinationPositionGenerator_SingleCrossover : public RecombinationPositionGenerator
{
    public:

    RecombinationPositionGenerator_SingleCrossover(const std::string& id)
    :   RecombinationPositionGenerator(id)
    {}

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "RecombinationPositionGenerator_SingleCrossover";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);

    private:
    std::vector<unsigned int> chromosome_lengths_;
};


//
// RecombinationPositionGenerator_Uniform
//

///
/// generates recombination events uniformly at random along the chromosome,
/// with a user-specified recombination rate for each chromosome
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// rate = \<float\> | 1 | optional; same rate (per chromosome) for all chromosomes
/// rates = \<float\> ... | none | optional; specific rates for all chromosomes, number of rates must match number of chromosomes
///
/// Note: 'rate' is the expected number of recombination positions (crossovers) for that chromosome per meiosis; 
///       the actual number of recombination positions ~ Poisson(rate).
///
/// Example: [tutorial_1_recombination_reporter.txt](../../examples/tutorial_1_recombination_reporter.txt)
///
/// \ingroup RecombinationPositionGenerators
///


class RecombinationPositionGenerator_Uniform : public RecombinationPositionGenerator
{
    public:

    struct ChromosomeInfo
    {
        size_t length;
        double rate;
        Random::DistributionPtr poisson;
        ChromosomeInfo(size_t _length, double _rate = 1.0);
    };

    RecombinationPositionGenerator_Uniform(const std::string& id, 
                                           std::vector<ChromosomeInfo> infos = std::vector<ChromosomeInfo>());

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "RecombinationPositionGenerator_Uniform";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);

    private:
    
    double common_rate_;
    typedef std::vector<ChromosomeInfo> ChromosomeInfos;
    ChromosomeInfos infos_;
};


//
// RecombinationPositionGenerator_RecombinationMap
//

///
/// generates recombination positions based on a genetic map
///
/// Note: genetic map files are expected to be in the HapMap genetic_map_* format
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// filename = \<genetic_map_filename\> | none | filename required for each chromosome pair
///
/// Example: [example_neutral_admixture.txt](../../examples/example_neutral_admixture.txt)
///
/// \ingroup RecombinationPositionGenerators
///


class RecombinationPositionGenerator_RecombinationMap : public RecombinationPositionGenerator
{ 
    public:

    RecombinationPositionGenerator_RecombinationMap(const std::string& id, 
        const std::vector<std::string>& filenames = std::vector<std::string>());

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "RecombinationPositionGenerator_RecombinationMap";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::vector<std::string> filenames_;
    std::vector< shared_ptr<RecombinationMap> > recombination_maps_;

    void read_files();
};


//
// RecombinationPositionGenerator_Composite
//

///
/// composite RecombinationPositionGenerator
///
/// parameter | default | notes
/// ----------|---------|-------------
/// chromosome:recombination_position_generator = \<chromosome\> \<id\> | none | optional
/// default_recombination_position_generator = \<id\> | RPG_Trivial | optional
///
/// Example: [example_sex_chromosome.txt](../../examples/example_sex_chromosome.txt)
///
/// \ingroup RecombinationPositionGenerators
///


class RecombinationPositionGenerator_Composite : public RecombinationPositionGenerator
{ 
    public:

    RecombinationPositionGenerator_Composite(const std::string& id);

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "RecombinationPositionGenerator_Composite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    RecombinationPositionGeneratorPtr default_rpg_;
    typedef std::map<size_t,RecombinationPositionGeneratorPtr> RPGMap;
    RPGMap rpg_map_;
};


#endif // _RECOMBINATIONPOSITIONGENERATORIMPLEMENTATION_HPP_

