//
// QuantitativeTraitImplementation.hpp
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


#ifndef _QUANTITATIVETRAITIMPLEMENTATION_HPP_
#define _QUANTITATIVETRAITIMPLEMENTATION_HPP_


#include "QuantitativeTrait.hpp"
#include "Random.hpp"


///
/// \defgroup QuantitativeTraits QuantitativeTraits
///
/// classes describing the genetic architecture of a quantitative trait; implementations
/// are responsible for the calculation of trait values from genotypes.
///


//
// QuantitativeTrait_PopulationComposite
//

///
/// per-population specification of quantitative traits
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_traits = \<id1\> [ \<id2\> ... ] | none | count must match population count
///
/// \ingroup QuantitativeTraits
///


class QuantitativeTrait_PopulationComposite : public QuantitativeTrait
{
    public:

    QuantitativeTrait_PopulationComposite(const std::string& id, 
                                          const QuantitativeTraitPtrs& qts = QuantitativeTraitPtrs()); 

    virtual void calculate_trait_values(const PopulationData& population_data) const;
    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_PopulationComposite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    QuantitativeTraitPtrs qts_;
};


//
// QuantitativeTrait_GenerationComposite
//

///
/// per-generation specification of quantitative traits
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// generation:quantitative_trait = \<int_generation_begin\> \<id_quantitative_trait\> | none | multiple allowed; generation 0 required
///
/// \ingroup QuantitativeTraits
///


class QuantitativeTrait_GenerationComposite : public QuantitativeTrait
{
    public:

    typedef std::map<size_t, QuantitativeTraitPtr> GenerationQTMap;

    QuantitativeTrait_GenerationComposite(const std::string& id, 
                                          const GenerationQTMap& qts = GenerationQTMap()); 

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_GenerationComposite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    GenerationQTMap qts_; // generation_index -> QT
};


//
// QuantitativeTrait_SingleLocusFitness 
//

///
/// represents fitness determined by a single locus
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// locus = \<id\> | none | required
/// w0 = \<float\> | none | required
/// w1 = \<float\> | none | required
/// w2 = \<float\> | none | required
///
/// Alternative parametrization:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus = \<id\> | none | required
/// additive_selection_coefficient = \<float\> | none | required; (w0,w1,w2) == (1,1+s,1+2s)
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup QuantitativeTraits
///

class QuantitativeTrait_SingleLocusFitness : public QuantitativeTrait
{
    public:

    QuantitativeTrait_SingleLocusFitness(const std::string& id, 
                                         Locus locus = Locus("id_dummy"), 
                                         std::vector<double> w = std::vector<double>(3));

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_SingleLocusFitness";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    Locus locus_;
    std::vector<double> w_; // relative fitnesses
};


//
// QTLEffect
//


struct QTLEffect
{
    Locus locus;
    std::vector<double> effects;

    QTLEffect();
    QTLEffect(const Locus& l, double e0, double e1, double e2);
    QTLEffect(const std::string& configuration, const Configurable::Registry& registry);
    std::string configuration() const;
};


typedef std::vector<QTLEffect> QTLEffects;


//
// QTLEffectGenerator
//

///
/// generates independent QTL effects based on specified distributions
/// 
/// Genotype effects:  0, a(1+k), 2a  
/// (where a is effect size, k is dominance -- Lynch & Walsh parametrization)
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus_list = \<id\> | none | required
/// effect_size_distribution = \<id\> | none | required
/// dominance_distribution = \<id\> | constant 0 | optional
///
/// \ingroup QuantitativeTraits
///

class QTLEffectGenerator : public Configurable
{
    public:

    QTLEffectGenerator(const std::string& id, 
                       LocusListPtr locus_list = LocusListPtr(),
                       Random::DistributionPtr effect_size_distribution = Random::DistributionPtr(),
                       Random::DistributionPtr dominance_distribution = Random::DistributionPtr(),
                       bool random_effect_sign = false);

    void generate_qtl_effects(QTLEffects& qtl_effects) const;

    // Configurable interface

    virtual std::string class_name() const {return "QTLEffectGenerator";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    LocusListPtr locus_list_;
    Random::DistributionPtr effect_size_distribution_;
    Random::DistributionPtr dominance_distribution_;
    bool random_effect_sign_;
};


typedef boost::shared_ptr<QTLEffectGenerator> QTLEffectGeneratorPtr;


//
// QuantitativeTrait_IndependentLoci
//

///
/// represents a trait determined by a multiple loci with independent effects
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// qtl = \<id_locus\> \<effect0\> \<effect1\> \<effect2\> | none | optional
/// qtl_effect_generator = \<id_qtl_effect_generator\> | none | optional
/// environmental_variance = \<double\> | 0 | optional
///
/// Example: [example_2_locus_selection.txt](../../examples/example_2_locus_selection.txt)
///
/// \ingroup QuantitativeTraits
///

class QuantitativeTrait_IndependentLoci : public QuantitativeTrait
{
    public:

    QuantitativeTrait_IndependentLoci(const std::string& id, 
                                      QTLEffects qtl_effects = QTLEffects(),
                                      double environmental_variance = 0);

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    const QTLEffects& qtl_effects() const {return qtl_effects_;}    

    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_IndependentLoci";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    QTLEffects qtl_effects_;
    size_t qtl_effects_specified_;
    double environmental_variance_;
    Random::DistributionPtr environment_effect_;
    std::vector<QTLEffectGeneratorPtr> qtl_effect_generators_;
};


//
// QuantitativeTrait_Expression
//

namespace mu {class Parser;}

///
/// represents a trait whose values are computed via a mathematical
/// expression from other trait values
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// expression = \<mathematical_expression\> | none | required
/// variable:quantitative_trait = \<variable\> \<id_quantitative_trait\> | none | optional
///
/// Example: [example_sex_chromosome.txt](../../examples/example_sex_chromosome.txt)
///
/// \ingroup QuantitativeTraits
///

class QuantitativeTrait_Expression : public QuantitativeTrait
{
    public:

    QuantitativeTrait_Expression(const std::string& id);

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_Expression";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    struct Assignment
    {
        std::string variable;
        std::string qtid;
    };

    typedef std::vector<Assignment> Assignments;
    Assignments assignments_;
    std::string expression_;

    boost::shared_ptr<mu::Parser> parser_;
};


//
// QuantitativeTrait_Alternator
//

///
/// assigns 0/1 alternately
/// 
/// Parameters: none
///
/// Example: TODO
///
/// \ingroup QuantitativeTraits
///

class QuantitativeTrait_Alternator : public QuantitativeTrait
{
    public:

    QuantitativeTrait_Alternator(const std::string& id) : QuantitativeTrait(id) {}

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    // Configurable interface

    virtual std::string class_name() const {return "QuantitativeTrait_Alternator";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}
};


#endif //  _QUANTITATIVETRAITIMPLEMENTATION_HPP_

