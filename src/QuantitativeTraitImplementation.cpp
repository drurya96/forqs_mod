//
// QuantitativeTraitImplementation.cpp
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


#include "QuantitativeTraitImplementation.hpp"
#include "Simulator.hpp"
#include "muparser/muParser.h"
#include <stdexcept>
#include <iostream>
#include <numeric>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


using namespace std;
namespace bfs = boost::filesystem;


// 
// QuantitativeTrait_PopulationComposite
//


QuantitativeTrait_PopulationComposite::QuantitativeTrait_PopulationComposite(const string& id, 
                                                                             const QuantitativeTraitPtrs& qts)
:   QuantitativeTrait(id), qts_(qts)
{}


void QuantitativeTrait_PopulationComposite::calculate_trait_values(
    const PopulationData& population_data, const Population& population) const
{
    if (population_data.population_index >= qts_.size())
        throw runtime_error("[QuantitativeTrait_PopulationComposite] Bad population_index.");

    QuantitativeTrait& qt = *qts_[population_data.population_index];
    TraitValueMap& trait_values = *population_data.trait_values;

    if (!trait_values.count(qt.object_id()))
        qt.calculate_trait_values(population_data, population);

    trait_values[object_id()] = trait_values[qt.object_id()];
}


Parameters QuantitativeTrait_PopulationComposite::parameters() const
{
    Parameters parameters;
    ostringstream qts_string;
    for (QuantitativeTraitPtrs::const_iterator it=qts_.begin(); it!=qts_.end(); ++it)
        qts_string << (*it)->object_id() << " ";
    parameters.insert_name_value("quantitative_traits", qts_string.str());
    return parameters;    
}


void QuantitativeTrait_PopulationComposite::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> qtids = parameters.value_vector<string>("quantitative_traits");
    for (vector<string>::const_iterator qtid=qtids.begin(); qtid!=qtids.end(); ++qtid)
        qts_.push_back(registry.get<QuantitativeTrait>(*qtid));
}


void QuantitativeTrait_PopulationComposite::initialize(const SimulatorConfig& config)
{
    if (qts_.size() != config.population_config_generator->population_count())
        throw runtime_error("[QuantitativeTrait_PopulationComposite] Quantitative trait count does not match population count.");
}


void QuantitativeTrait_PopulationComposite::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (QuantitativeTraitPtrs::const_iterator it=qts_.begin(); it!=qts_.end(); ++it)
    {
        (*it)->write_configuration(os, ids_written);
        ids_written.insert((*it)->object_id());
    }
}


// 
// QuantitativeTrait_GenerationComposite
//


QuantitativeTrait_GenerationComposite::QuantitativeTrait_GenerationComposite(const string& id, 
                                                                             const GenerationQTMap& qts)
:   QuantitativeTrait(id), qts_(qts)
{}


void QuantitativeTrait_GenerationComposite::calculate_trait_values(const PopulationData& population_data, const Population& population) const
{
    if (qts_.empty())
        throw runtime_error("[QuantitativeTrait_GenerationComposite] Initialization error: no quantitative traits."
                            "Check parameter name: generation:quantitative_trait");

    if (!qts_.count(0))
        throw runtime_error("[QuantitativeTrait_GenerationComposite] Initialization error: generation 0 not specified.");

    GenerationQTMap::const_iterator it = qts_.upper_bound(population_data.generation_index);
    if (it != qts_.begin()) --it;

    const QuantitativeTrait& qt = *it->second;
    
    TraitValueMap& trait_values = *population_data.trait_values;

    if (!trait_values.count(qt.object_id()))
        qt.calculate_trait_values(population_data, population);

    trait_values[object_id()] = trait_values[qt.object_id()];
}


Parameters QuantitativeTrait_GenerationComposite::parameters() const
{
    Parameters parameters;
    for (GenerationQTMap::const_iterator it=qts_.begin(); it!=qts_.end(); ++it)
        parameters.insert_name_value("generation:quantitative_trait", boost::lexical_cast<string>(it->first) + " " + it->second->object_id());
    return parameters;    
}


void QuantitativeTrait_GenerationComposite::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> parameter_values = parameters.values<string>("generation:quantitative_trait");    

    for (vector<string>::const_iterator it=parameter_values.begin(); it!=parameter_values.end(); ++it)
    {
        istringstream iss(*it);
        size_t generation_index;
        string qtid;
        iss >> generation_index >> qtid;
        qts_[generation_index] = registry.get<QuantitativeTrait>(qtid);
    }
}


void QuantitativeTrait_GenerationComposite::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (GenerationQTMap::const_iterator it=qts_.begin(); it!=qts_.end(); ++it)
    {
        it->second->write_configuration(os, ids_written);
        ids_written.insert(it->second->object_id());
    }
}


// 
// QuantitativeTrait_SingleLocusFitness
//


QuantitativeTrait_SingleLocusFitness::QuantitativeTrait_SingleLocusFitness(const string& id, Locus locus, vector<double> w)
:   QuantitativeTrait(id), locus_(locus), w_(w)
{
    loci_.insert(locus);
}


void QuantitativeTrait_SingleLocusFitness::calculate_trait_values(const PopulationData& population_data, const Population& population) const
{
    const GenotypeDataPtr& g = population_data.genotypes->get(locus_);

    DataVectorPtr fitnesses(new DataVector(g->size()));

    // transform {0, 1, 2} -> {w[0], w[1], w[2]}
    DataVector::iterator jt = fitnesses->begin();
    for (GenotypeData::const_iterator it=g->begin(); it!=g->end(); ++it, ++jt)
        *jt = w_[genotype_sum(*it)];

    (*population_data.trait_values)[object_id()] = fitnesses;
}


Parameters QuantitativeTrait_SingleLocusFitness::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("locus", locus_.object_id());
    parameters.insert_name_value("w0", w_[0]);
    parameters.insert_name_value("w1", w_[1]);
    parameters.insert_name_value("w2", w_[2]);
    return parameters;
}


void QuantitativeTrait_SingleLocusFitness::configure(const Parameters& parameters, 
                                                     const Registry& registry)
{
    loci_.clear();
    locus_ = *registry.get<Locus>(parameters.value<string>("locus"));
    loci_.insert(locus_);

    if (parameters.count("additive_selection_coefficient"))
    {
        double s = parameters.value<double>("additive_selection_coefficient");
        w_[0] = 1;
        w_[1] = 1 + s;
        w_[2] = 1 + 2*s;
    }
    else
    {
        w_[0] = parameters.value<double>("w0");
        w_[1] = parameters.value<double>("w1");
        w_[2] = parameters.value<double>("w2");
    }
}


void QuantitativeTrait_SingleLocusFitness::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    locus_.write_configuration(os, ids_written);
}


// 
// QTLEffect
//


QTLEffect::QTLEffect()
:   locus("id_dummy"), effects(3) 
{}


QTLEffect::QTLEffect(const Locus& l, double e0, double e1, double e2)
:   locus(l), effects(3)
{
    effects[0] = e0;
    effects[1] = e1;
    effects[2] = e2;
}


QTLEffect::QTLEffect(const string& configuration, const Configurable::Registry& registry)
:   locus("id_dummy"), effects(3)
{
    string locus_id;
    istringstream iss(configuration);
    iss >> locus_id >> effects[0] >> effects[1] >> effects[2];
    locus = *registry.get<Locus>(locus_id);
}


string QTLEffect::configuration() const
{
    ostringstream result;
    result << locus.object_id() << " " << effects[0] << " " << effects[1] << " " << effects[2];
    return result.str();
}


//
// QTLEffectGenerator
//


QTLEffectGenerator::QTLEffectGenerator(const string& id, 
                                       LocusListPtr locus_list,
                                       Random::DistributionPtr effect_size_distribution,
                                       Random::DistributionPtr dominance_distribution,
                                       bool random_effect_sign)
:   Configurable(id), 
    locus_list_(locus_list),
    effect_size_distribution_(effect_size_distribution),
    dominance_distribution_(dominance_distribution),
    random_effect_sign_(random_effect_sign)
{}


void QTLEffectGenerator::generate_qtl_effects(QTLEffects& qtl_effects) const
{
    if (!locus_list_.get())
        throw runtime_error("[QTLEffectGenerator] Null LocusList.");

    if (!effect_size_distribution_.get())
        throw runtime_error("[QTLEffectGenerator] Null effect size distribution.");

    if (!dominance_distribution_.get())
        throw runtime_error("[QTLEffectGenerator] Null dominance distribution.");

    for (LocusList::const_iterator locus=locus_list_->begin(); locus!=locus_list_->end(); ++locus)
    {
        double a = effect_size_distribution_->random_value();
        if (random_effect_sign_ && Random::bernoulli()) a = -a;
        double k = dominance_distribution_->random_value();
        qtl_effects.push_back(QTLEffect(*locus, 0, a*(1+k), 2*a));
    }
}


Parameters QTLEffectGenerator::parameters() const
{
    Parameters parameters;

    if (locus_list_.get())
        parameters.insert_name_value("locus_list", locus_list_->object_id());

    if (effect_size_distribution_.get())
        parameters.insert_name_value("effect_size_distribution", effect_size_distribution_->object_id());

    if (random_effect_sign_)
        parameters.insert_name_value("random_effect_sign", random_effect_sign_);

    if (dominance_distribution_.get())
        parameters.insert_name_value("dominance_distribution", dominance_distribution_->object_id());

    return parameters;
}


void QTLEffectGenerator::configure(const Parameters& parameters, const Registry& registry)
{
    locus_list_ = 
        registry.get<LocusList>(parameters.value<string>("locus_list"));

    effect_size_distribution_ = 
        registry.get<Random::Distribution>(parameters.value<string>("effect_size_distribution"));

    random_effect_sign_ = parameters.value<bool>("random_effect_sign", false);

    if (parameters.count("dominance_distribution"))
        dominance_distribution_ = registry.get<Random::Distribution>(parameters.value<string>("dominance_distribution"));
    else
        dominance_distribution_ = Random::create_constant_distribution("id_dummy_zero", 0.);
}


void QTLEffectGenerator::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    if (locus_list_.get()) locus_list_->write_configuration(os, ids_written);
    if (effect_size_distribution_.get()) effect_size_distribution_->write_configuration(os, ids_written);
    if (dominance_distribution_.get()) dominance_distribution_->write_configuration(os, ids_written);
}


// 
// QuantitativeTrait_IndependentLoci
//


QuantitativeTrait_IndependentLoci::QuantitativeTrait_IndependentLoci(const std::string& id, 
                                                                     QTLEffects qtl_effects,
                                                                     double environmental_variance)
:   QuantitativeTrait(id), 
    qtl_effects_(qtl_effects), 
    qtl_effects_specified_(qtl_effects.size()),
    environmental_variance_(environmental_variance)
{}


void QuantitativeTrait_IndependentLoci::calculate_trait_values(const PopulationData& population_data, const Population& population) const
{
    DataVectorPtr trait_values(new DataVector(population_data.population_size));
    (*population_data.trait_values)[object_id()] = trait_values;

    // sanity check: make sure we have the necessary genotype data

    for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.end(); ++it)
    {
        if (!population_data.genotypes->count(it->locus))
        {
            ostringstream oss;
            oss << "[QuantitativeTrait_IndependentLoci] Invalid genotype map: locus "
                << it->locus << " not found.";
            throw runtime_error(oss.str().c_str());
        }
    }

    // calculate trait values

    for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.end(); ++it)
    {
        const GenotypeData& g = *population_data.genotypes->get(it->locus);

        if (g.size() != population_data.population_size)
            throw runtime_error("[QuantitativeTrait_IndependentLoci] Genotype vector sizes don't match.");

        // transform genotype (gt) to effect size:  
        //     {0, 1, 2} -> {effects[0], effects[1], effects[2]}
        // and update trait value (tv) by the effect size

        DataVector::iterator tv = trait_values->begin();
        for (GenotypeData::const_iterator gt=g.begin(); gt!=g.end(); ++gt, ++tv)
            *tv += it->effects[genotype_sum(*gt)];
    }

    if (environment_effect_.get())
        for (DataVector::iterator tv=trait_values->begin(); tv!=trait_values->end(); ++tv)
            *tv += environment_effect_->random_value();
}


Parameters QuantitativeTrait_IndependentLoci::parameters() const
{
    Parameters parameters;

    for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.begin() + qtl_effects_specified_; ++it)
        parameters.insert_name_value("qtl", it->configuration());

    for (vector<QTLEffectGeneratorPtr>::const_iterator gen=qtl_effect_generators_.begin(); 
         gen!=qtl_effect_generators_.end(); ++gen)
        parameters.insert_name_value("qtl_effect_generator", (*gen)->object_id());

    parameters.insert_name_value("environmental_variance", environmental_variance_);

    return parameters;
}


void QuantitativeTrait_IndependentLoci::configure(const Parameters& parameters, 
                                                  const Registry& registry)
{
    qtl_effects_.clear();

    // handle specified qtls and environmental variance

    vector<string> qtl_configurations = parameters.values<string>("qtl");

    for (vector<string>::const_iterator configuration=qtl_configurations.begin();
            configuration!=qtl_configurations.end(); ++configuration)
        qtl_effects_.push_back(QTLEffect(*configuration, registry));

    qtl_effects_specified_ = qtl_configurations.size();

    environmental_variance_ = parameters.value<double>("environmental_variance", 0.);

    // save references to QTLEffectGenerators

    vector<string> qtl_effect_generator_ids = parameters.values<string>("qtl_effect_generator");

    for (vector<string>::const_iterator id=qtl_effect_generator_ids.begin(); id!=qtl_effect_generator_ids.end(); ++id)
        qtl_effect_generators_.push_back(registry.get<QTLEffectGenerator>(*id));
}


void QuantitativeTrait_IndependentLoci::initialize(const SimulatorConfig& config)
{
    // generate random QTL effects

    for (vector<QTLEffectGeneratorPtr>::const_iterator gen=qtl_effect_generators_.begin(); 
         gen!=qtl_effect_generators_.end(); ++gen)
        (*gen)->generate_qtl_effects(qtl_effects_);

    // other initialization

    loci_.clear();

    for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.end(); ++it)
        loci_.insert(it->locus);

    if (environmental_variance_ > 0) 
        environment_effect_ = Random::create_normal_distribution(
            "QT_IndependentLoci_environmental_effect_", 0, environmental_variance_);
    else
        environment_effect_ = Random::DistributionPtr();

    // write output file with QTL info

    if (!config.output_directory.empty())
    {
        bfs::path output_directory = config.output_directory;
        bfs::path filename = output_directory / "qt_independent_loci.qtl.txt";

        int index = 1;
        while (bfs::exists(filename))
        {
            ++index;
            ostringstream oss;
            oss << "qt_independent_loci.qtl" << index << ".txt";
            filename = output_directory / oss.str();
        }

        bfs::ofstream os(filename);
        if (!os)
            throw runtime_error(("[QuantitativeTrait_IndependentLoci] Unable to open file " + filename.string()).c_str());

        os << "# id chromosome position e0 e1 e2\n";
        for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.end(); ++it)
            os << it->locus.object_id() << " "
               << it->locus.chromosome_pair_index + 1 << " "
               << it->locus.position << " "
               << it->effects[0] << " "
               << it->effects[1] << " "
               << it->effects[2] << endl;
    }
}


void QuantitativeTrait_IndependentLoci::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (QTLEffects::const_iterator it=qtl_effects_.begin(); it!=qtl_effects_.begin() + qtl_effects_specified_; ++it)
        it->locus.write_configuration(os, ids_written);

    for (vector<QTLEffectGeneratorPtr>::const_iterator gen=qtl_effect_generators_.begin(); 
         gen!=qtl_effect_generators_.end(); ++gen)
        (*gen)->write_configuration(os, ids_written);
}


//
// QuantitativeTrait_Expression
//


QuantitativeTrait_Expression::QuantitativeTrait_Expression(const string& id)
:   QuantitativeTrait(id)
{}


void QuantitativeTrait_Expression::calculate_trait_values(const PopulationData& population_data, const Population& population) const
{
    try
    {
        DataVectorPtr result(new DataVector(population_data.population_size));
        (*population_data.trait_values)[object_id()] = result;

        if (result->empty()) return;

        for (Assignments::const_iterator assignment = assignments_.begin();
            assignment != assignments_.end(); ++assignment)
        {
            DataVector& data = *population_data.trait_values->get(assignment->qtid);
            if (data.empty())
                throw runtime_error("[QuantitativeTrait_Expression::calculate_trait_values()] Empty data.");
            parser_->DefineVar(assignment->variable, &data[0]);
        }

        parser_->Eval(&(*result)[0], population_data.population_size);
    }
    catch (mu::ParserError& e)
    {
        ostringstream message;
        message << "[QuantitativeTrait_Expression::calculate_trait_values()] Exception thrown from muparser library:\n"
                << e.GetMsg();
        throw runtime_error(message.str());
    }
}


Parameters QuantitativeTrait_Expression::parameters() const
{
    Parameters parameters;

    for (Assignments::const_iterator assignment = assignments_.begin();
        assignment != assignments_.end(); ++assignment)
        parameters.insert_name_value("variable:quantitative_trait",
            assignment->variable + " " + assignment->qtid);

    parameters.insert_name_value("expression", expression_);

    return parameters;
}


void QuantitativeTrait_Expression::configure(const Parameters& parameters, const Registry& registry)
{
    assignments_.clear();

    vector<string> assignment_strings = parameters.values<string>("variable:quantitative_trait");
    for (vector<string>::const_iterator assignment_string = assignment_strings.begin();
        assignment_string != assignment_strings.end(); ++assignment_string)
    {
        istringstream iss(*assignment_string);
        Assignment assignment;
        iss >> assignment.variable >> assignment.qtid;
        assignments_.push_back(assignment);
    }

    expression_ = parameters.value<string>("expression");

    try
    {
        parser_ = boost::shared_ptr<mu::Parser>(new mu::Parser);
        parser_->SetExpr(expression_);
    }
    catch (mu::ParserError& e)
    {
        ostringstream message;
        message << "[QuantitativeTrait_Expression::configure()] Exception thrown from muparser library:\n"
                << e.GetMsg();
        throw runtime_error(message.str());
    }
}


//
// QuantitativeTrait_Alternator
//

void QuantitativeTrait_Alternator::calculate_trait_values(const PopulationData& population_data, const Population& population) const
{
    DataVectorPtr result(new DataVector(population_data.population_size));
    (*population_data.trait_values)[object_id()] = result;

    for (size_t i=0; i<result->size(); ++i)
        (*result)[i] = i%2;
}


