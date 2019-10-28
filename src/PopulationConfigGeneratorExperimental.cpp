//
// PopulationConfigGeneratorExperimental.cpp
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


#include "PopulationConfigGeneratorExperimental.hpp"
#include "QuantitativeTraitImplementation.hpp"
#include "FitnessFunctionImplementation.hpp"
#include "VariantIndicatorImplementation.hpp"
#include "ReporterImplementation.hpp"
#include "Simulator.hpp"
#include "boost/lambda/casts.hpp"
#include <fstream>
#include <stdexcept>
#include <numeric>


using namespace std;


//
// PopulationConfigGenerator_TurnerExperiment
//


PopulationConfigGenerator_TurnerExperiment::PopulationConfigGenerator_TurnerExperiment(const string& id)
:   PopulationConfigGenerator(id), 
    population_size_founders_(0),
    population_size_neutral_(0),
    population_size_replicate_(0),
    generation_count_neutral_(0),
    generation_count_selected_(0),
    proportion_selected_(0), 
    report_hidden_(false)
{
    // hard-coded Drosophila chromsome lengths
    chromosome_lengths_.clear();
    chromosome_lengths_.push_back(2.25e7);  // X
    chromosome_lengths_.push_back(4.62e7);  // 2L: 0-23.1mb, 2R: 25-46.2mb
    chromosome_lengths_.push_back(5.3e7);   // 3L: 0-24.6mb, 3R: 25-52.9mb
    chromosome_pair_count_ = chromosome_lengths_.size();
}


Population::Configs PopulationConfigGenerator_TurnerExperiment::population_configs(
    size_t generation_index, const PopulationDataPtrs& population_datas) const
{
    // note: trait names follow the "user" convention, i.e. population number is 1-based;
    //       however, MatingDistribution::Entry uses 0-based population indices

    const size_t pop1 = 0;
    const size_t pop2 = 1;
    const size_t pop3 = 2;
    const size_t pop4 = 3;

    if (generation_index == 0)
    {
        Population::Configs popconfigs(1);
        Population::Config& popconfig = popconfigs.front();
        popconfig.chromosome_pair_count = chromosome_pair_count_;
        popconfig.population_size = population_size_founders_;
        return popconfigs;
    }
    else if (generation_index <= generation_count_neutral_)
    {
        // random mating

        Population::Configs popconfigs(1);
        Population::Config& popconfig = popconfigs.front();
        popconfig.chromosome_pair_count = chromosome_pair_count_;
        popconfig.population_size = population_size_neutral_;
        popconfig.mating_distribution.push_back(
            MatingDistribution::Entry(1, pop1, "_female", pop1, "_male"));
        return popconfigs;
    }
    else if (generation_index == generation_count_neutral_ + 1)
    {
        // first round of selection: split into 4 populations

        Population::Configs popconfigs(4);

        for (Population::Configs::iterator popconfig=popconfigs.begin(); 
             popconfig!=popconfigs.end(); ++popconfig)
        {
            popconfig->chromosome_pair_count = chromosome_pair_count_;
            popconfig->population_size = population_size_replicate_;
        }

        popconfigs[0].mating_distribution.push_back(
            MatingDistribution::Entry(1, pop1, "_female", pop1, "_selected_low_1_even"));
        popconfigs[1].mating_distribution.push_back(
            MatingDistribution::Entry(1, pop1, "_female", pop1, "_selected_low_1_odd"));
        popconfigs[2].mating_distribution.push_back(
            MatingDistribution::Entry(1, pop1, "_female", pop1, "_selected_high_1_even"));
        popconfigs[3].mating_distribution.push_back(
            MatingDistribution::Entry(1, pop1, "_female", pop1, "_selected_high_1_odd"));

        return popconfigs;
    }
    else if (generation_index <= generation_count_)
    {
        // subsequent rounds of selection with migration

        if (population_datas.size() != 4)
            throw runtime_error("[PopulationConfigGenerator_TurnerExperiment] Bad population_datas.");

        Population::Configs popconfigs(4);

        for (Population::Configs::iterator popconfig=popconfigs.begin(); 
             popconfig!=popconfigs.end(); ++popconfig)
        {
            popconfig->chromosome_pair_count = chromosome_pair_count_;
            popconfig->population_size = population_size_replicate_;
        }

        // migration rates are given by absolute counts of selected individuals

        // population 1

        vector<size_t> selected_counts(3);

        selected_counts[0] = population_datas[pop1]->trait_values->get("_selected_low_1")->sum();
        selected_counts[1] = population_datas[pop3]->trait_values->get("_selected_low_1_even")->sum();
        selected_counts[2] = population_datas[pop4]->trait_values->get("_selected_low_1_even")->sum();

        if (selected_counts[0]) popconfigs[0].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[0], pop1, "_female", pop1, "_selected_low_1"));
        if (selected_counts[1]) popconfigs[0].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[1], pop1, "_female", pop3, "_selected_low_1_even"));
        if (selected_counts[2]) popconfigs[0].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[2], pop1, "_female", pop4, "_selected_low_1_even"));

        // population 2

        selected_counts[0] = population_datas[pop2]->trait_values->get("_selected_low_2")->sum();
        selected_counts[1] = population_datas[pop3]->trait_values->get("_selected_low_2_odd")->sum();
        selected_counts[2] = population_datas[pop4]->trait_values->get("_selected_low_2_odd")->sum();

        if (selected_counts[0]) popconfigs[1].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[0], pop2, "_female", pop2, "_selected_low_2"));
        if (selected_counts[1]) popconfigs[1].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[1], pop2, "_female", pop3, "_selected_low_2_odd"));
        if (selected_counts[2]) popconfigs[1].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[2], pop2, "_female", pop4, "_selected_low_2_odd"));

        // population 3

        selected_counts[0] = population_datas[pop3]->trait_values->get("_selected_high_3")->sum();
        selected_counts[1] = population_datas[pop1]->trait_values->get("_selected_high_3_even")->sum();
        selected_counts[2] = population_datas[pop2]->trait_values->get("_selected_high_3_even")->sum();

        if (selected_counts[0]) popconfigs[2].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[0], pop3, "_female", pop3, "_selected_high_3"));
        if (selected_counts[1]) popconfigs[2].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[1], pop3, "_female", pop1, "_selected_high_3_even"));
        if (selected_counts[2]) popconfigs[2].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[2], pop3, "_female", pop2, "_selected_high_3_even"));

        // population 4

        selected_counts[0] = population_datas[pop4]->trait_values->get("_selected_high_4")->sum();
        selected_counts[1] = population_datas[pop1]->trait_values->get("_selected_high_4_odd")->sum();
        selected_counts[2] = population_datas[pop2]->trait_values->get("_selected_high_4_odd")->sum();

        if (selected_counts[0]) popconfigs[3].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[0], pop4, "_female", pop4, "_selected_high_4"));
        if (selected_counts[1]) popconfigs[3].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[1], pop4, "_female", pop1, "_selected_high_4_odd"));
        if (selected_counts[2]) popconfigs[3].mating_distribution.push_back(
            MatingDistribution::Entry(selected_counts[2], pop4, "_female", pop2, "_selected_high_4_odd"));

        return popconfigs;
    }
    else
    {
        throw runtime_error("[PopulationConfigGenerator_TurnerExperiment] Invalid generation index.");
    }
}


Parameters PopulationConfigGenerator_TurnerExperiment::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("population_size_founders", population_size_founders_);
    parameters.insert_name_value("population_size_neutral", population_size_neutral_);
    parameters.insert_name_value("population_size_replicate", population_size_replicate_);
    parameters.insert_name_value("generation_count_neutral", generation_count_neutral_);
    parameters.insert_name_value("generation_count_selected", generation_count_selected_);
    parameters.insert_name_value("quantitative_trait", qtid_);
    parameters.insert_name_value("proportion_selected", proportion_selected_);
    parameters.insert_name_value("report_hidden", report_hidden_);
    return parameters;
}


void PopulationConfigGenerator_TurnerExperiment::configure(const Parameters& parameters, const Registry& registry)
{
    population_size_founders_ = parameters.value<size_t>("population_size_founders");
    population_size_neutral_ = parameters.value<size_t>("population_size_neutral");
    population_size_replicate_ = parameters.value<size_t>("population_size_replicate");
    generation_count_neutral_ = parameters.value<size_t>("generation_count_neutral");
    generation_count_selected_ = parameters.value<size_t>("generation_count_selected");
    qtid_ = parameters.value<string>("quantitative_trait");
    proportion_selected_ = parameters.value<double>("proportion_selected");
    report_hidden_ = parameters.value<bool>("report_hidden", false);

    generation_count_ = generation_count_neutral_ + generation_count_selected_;
}


void PopulationConfigGenerator_TurnerExperiment::initialize(const SimulatorConfig& config)
{
    SimulatorConfig& simconfig = const_cast<SimulatorConfig&>(config);

    Configurable::Registry registry;

    // gender

    LocusPtr locus_Y(new Locus("_locus_Y", 0, 0)); // 0-based chromosome index, position
    registry[locus_Y->object_id()] = locus_Y;

    /*
    LocusListPtr locus_list_Y(new LocusList("_locus_list_Y"));
    locus_list_Y->push_back(*locus_Y);
    registry[locus_list_Y->object_id()] = locus_list_Y;
    */

    QuantitativeTraitPtr male(new QuantitativeTrait_IndependentLoci("_male"));
    Parameters parameters_male;
    parameters_male.insert_name_value("qtl", "_locus_Y 0 1 2");
    male->configure(parameters_male, registry);
    male->initialize(simconfig);
    simconfig.quantitative_traits.push_back(male);

    QuantitativeTraitPtr female(new QuantitativeTrait_Expression("_female"));
    Parameters parameters_female;
    parameters_female.insert_name_value("variable:quantitative_trait", "y _male");
    parameters_female.insert_name_value("expression", "1-y");
    female->configure(parameters_female, registry);
    simconfig.quantitative_traits.push_back(female);

    /*
    // automatic male/female specification

    VariantIndicatorPtr vi_Y(new VariantIndicator_IDRange("_vi_Y"));
    registry[vi_Y->object_id()] = vi_Y;
    size_t haplotype_count = population_size_founders_ * 2;
    ostringstream vi_Y_value;
    vi_Y_value << "_locus_Y 1 " << haplotype_count << " 4 1"; // alternate male/female
    Parameters parameters_vi_Y;
    parameters_vi_Y.insert_name_value("locus:start:count:step:value", vi_Y_value.str());
    vi_Y->configure(parameters_vi_Y, registry);

    VariantIndicatorPtr vi_composite(new VariantIndicator_Composite("_vi_composite"));
    registry[simconfig.variant_indicator->object_id()] = simconfig.variant_indicator;
    ostringstream vi_composite_value;
    vi_composite_value << simconfig.variant_indicator->object_id() << " " << vi_Y->object_id();
    Parameters parameters_vi_composite;
    parameters_vi_composite.insert_name_value("variant_indicators", vi_composite_value.str());
    vi_composite->configure(parameters_vi_composite, registry);
    simconfig.variant_indicator = vi_composite;
    */

    // selection

    QuantitativeTraitPtr trait_male(new QuantitativeTrait_Expression("_trait_male"));
    Parameters parameters_trait_male;
    parameters_trait_male.insert_name_value("variable:quantitative_trait", "trait " + qtid_);
    parameters_trait_male.insert_name_value("variable:quantitative_trait", "male _male");
    parameters_trait_male.insert_name_value("expression", "trait * male");
    trait_male->configure(parameters_trait_male, registry);
    simconfig.quantitative_traits.push_back(trait_male);

    Parameters parameters_selected_common;
    parameters_selected_common.insert_name_value("quantitative_trait", "_trait_male");
    parameters_selected_common.insert_name_value("proportion_selected", proportion_selected_);
    parameters_selected_common.insert_name_value("ignore_zero", 1);

    QuantitativeTraitPtr selected_low_1(new FitnessFunction_TruncationSelection("_selected_low_1"));
    Parameters parameters_selected_low_1 = parameters_selected_common;
    parameters_selected_low_1.insert_name_value("single_threshold_population", 1);
    parameters_selected_low_1.insert_name_value("lower_tail", 1);
    selected_low_1->configure(parameters_selected_low_1, registry);
    simconfig.quantitative_traits.push_back(selected_low_1);

    QuantitativeTraitPtr selected_high_1(new FitnessFunction_TruncationSelection("_selected_high_1"));
    Parameters parameters_selected_high_1 = parameters_selected_common;
    parameters_selected_high_1.insert_name_value("single_threshold_population", 1);
    selected_high_1->configure(parameters_selected_high_1, registry);
    simconfig.quantitative_traits.push_back(selected_high_1);

    QuantitativeTraitPtr selected_low_2(new FitnessFunction_TruncationSelection("_selected_low_2"));
    Parameters parameters_selected_low_2 = parameters_selected_common;
    parameters_selected_low_2.insert_name_value("single_threshold_population", 2);
    parameters_selected_low_2.insert_name_value("lower_tail", 1);
    selected_low_2->configure(parameters_selected_low_2, registry);
    simconfig.quantitative_traits.push_back(selected_low_2);

    QuantitativeTraitPtr selected_high_3(new FitnessFunction_TruncationSelection("_selected_high_3"));
    Parameters parameters_selected_high_3 = parameters_selected_common;
    parameters_selected_high_3.insert_name_value("single_threshold_population", 3);
    selected_high_3->configure(parameters_selected_high_3, registry);
    simconfig.quantitative_traits.push_back(selected_high_3);

    QuantitativeTraitPtr selected_high_4(new FitnessFunction_TruncationSelection("_selected_high_4"));
    Parameters parameters_selected_high_4 = parameters_selected_common;
    parameters_selected_high_4.insert_name_value("single_threshold_population", 4);
    selected_high_4->configure(parameters_selected_high_4, registry);
    simconfig.quantitative_traits.push_back(selected_high_4);

    // migration

    QuantitativeTraitPtr even(new QuantitativeTrait_Alternator("_even"));
    simconfig.quantitative_traits.push_back(even);

    QuantitativeTraitPtr odd(new QuantitativeTrait_Expression("_odd"));
    Parameters parameters_odd;
    parameters_odd.insert_name_value("variable:quantitative_trait", "even _even");
    parameters_odd.insert_name_value("expression", "1 - even");
    odd->configure(parameters_odd, registry);
    simconfig.quantitative_traits.push_back(odd);

    // split; migrants: 3,4 -> 1
    QuantitativeTraitPtr selected_low_1_even(new QuantitativeTrait_Expression("_selected_low_1_even"));
    Parameters parameters_selected_low_1_even;
    parameters_selected_low_1_even.insert_name_value("variable:quantitative_trait", "selected_low_1 _selected_low_1");
    parameters_selected_low_1_even.insert_name_value("variable:quantitative_trait", "even _even");
    parameters_selected_low_1_even.insert_name_value("expression", "selected_low_1 * even");
    selected_low_1_even->configure(parameters_selected_low_1_even, registry);
    simconfig.quantitative_traits.push_back(selected_low_1_even);

    // split
    QuantitativeTraitPtr selected_low_1_odd(new QuantitativeTrait_Expression("_selected_low_1_odd"));
    Parameters parameters_selected_low_1_odd;
    parameters_selected_low_1_odd.insert_name_value("variable:quantitative_trait", "selected_low_1 _selected_low_1");
    parameters_selected_low_1_odd.insert_name_value("variable:quantitative_trait", "odd _odd");
    parameters_selected_low_1_odd.insert_name_value("expression", "selected_low_1 * odd");
    selected_low_1_odd->configure(parameters_selected_low_1_odd, registry);
    simconfig.quantitative_traits.push_back(selected_low_1_odd);

    // split
    QuantitativeTraitPtr selected_high_1_even(new QuantitativeTrait_Expression("_selected_high_1_even"));
    Parameters parameters_selected_high_1_even;
    parameters_selected_high_1_even.insert_name_value("variable:quantitative_trait", "selected_high_1 _selected_high_1");
    parameters_selected_high_1_even.insert_name_value("variable:quantitative_trait", "even _even");
    parameters_selected_high_1_even.insert_name_value("expression", "selected_high_1 * even");
    selected_high_1_even->configure(parameters_selected_high_1_even, registry);
    simconfig.quantitative_traits.push_back(selected_high_1_even);

    // split
    QuantitativeTraitPtr selected_high_1_odd(new QuantitativeTrait_Expression("_selected_high_1_odd"));
    Parameters parameters_selected_high_1_odd;
    parameters_selected_high_1_odd.insert_name_value("variable:quantitative_trait", "selected_high_1 _selected_high_1");
    parameters_selected_high_1_odd.insert_name_value("variable:quantitative_trait", "odd _odd");
    parameters_selected_high_1_odd.insert_name_value("expression", "selected_high_1 * odd");
    selected_high_1_odd->configure(parameters_selected_high_1_odd, registry);
    simconfig.quantitative_traits.push_back(selected_high_1_odd);

    // migrants: 3,4 -> 2
    QuantitativeTraitPtr selected_low_2_odd(new QuantitativeTrait_Expression("_selected_low_2_odd"));
    Parameters parameters_selected_low_2_odd;
    parameters_selected_low_2_odd.insert_name_value("variable:quantitative_trait", "selected_low_2 _selected_low_2");
    parameters_selected_low_2_odd.insert_name_value("variable:quantitative_trait", "odd _odd");
    parameters_selected_low_2_odd.insert_name_value("expression", "selected_low_2 * odd");
    selected_low_2_odd->configure(parameters_selected_low_2_odd, registry);
    simconfig.quantitative_traits.push_back(selected_low_2_odd);

    // migrants: 1,2 -> 3
    QuantitativeTraitPtr selected_high_3_even(new QuantitativeTrait_Expression("_selected_high_3_even"));
    Parameters parameters_selected_high_3_even;
    parameters_selected_high_3_even.insert_name_value("variable:quantitative_trait", "selected_high_3 _selected_high_3");
    parameters_selected_high_3_even.insert_name_value("variable:quantitative_trait", "even _even");
    parameters_selected_high_3_even.insert_name_value("expression", "selected_high_3 * even");
    selected_high_3_even->configure(parameters_selected_high_3_even, registry);
    simconfig.quantitative_traits.push_back(selected_high_3_even);

    // migrants: 1,2 -> 4
    QuantitativeTraitPtr selected_high_4_odd(new QuantitativeTrait_Expression("_selected_high_4_odd"));
    Parameters parameters_selected_high_4_odd;
    parameters_selected_high_4_odd.insert_name_value("variable:quantitative_trait", "selected_high_4 _selected_high_4");
    parameters_selected_high_4_odd.insert_name_value("variable:quantitative_trait", "odd _odd");
    parameters_selected_high_4_odd.insert_name_value("expression", "selected_high_4 * odd");
    selected_high_4_odd->configure(parameters_selected_high_4_odd, registry);
    simconfig.quantitative_traits.push_back(selected_high_4_odd);

    // reporter

    if (report_hidden_)
    {
        ReporterPtr reporter_trait_values(new Reporter_TraitValues("_reporter_trait_values"));
        Parameters parameters_reporter_trait_values;
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_male _female");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", qtid_);
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_trait_male");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_selected_low_1 _selected_high_1");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_selected_low_2 _selected_high_3 _selected_high_4");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_selected_low_1_even _selected_low_1_odd");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_selected_high_1_even _selected_high_1_odd");
        parameters_reporter_trait_values.insert_name_value("quantitative_traits", "_selected_low_2_odd _selected_high_3_even _selected_high_4_odd");
        parameters_reporter_trait_values.insert_name_value("ignore_zero_values", qtid_ + " _trait_male");
        parameters_reporter_trait_values.insert_name_value("write_full", 1);
        parameters_reporter_trait_values.insert_name_value("filetag", "hidden");
        reporter_trait_values->configure(parameters_reporter_trait_values, registry);
        reporter_trait_values->initialize(simconfig);
        simconfig.reporters.push_back(reporter_trait_values);
    }
}



