//
// QuantitativeTraitTestImplementation.cpp
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
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_Configurable_QuantitativeTrait_SingleLocusFitness()
{
    if (os_) *os_ << "test_Configurable_QuantitativeTrait_SingleLocusFitness()\n";

    string locus_id = "my_locus";
    LocusPtr locus(new Locus(locus_id, 0, 123456));

    Configurable::Registry registry;
    registry[locus_id] = locus;

    Parameters parameters_in;
    parameters_in.insert_name_value("locus", locus_id);
    parameters_in.insert_name_value("w0", 1.0);
    parameters_in.insert_name_value("w1", 1.1);
    parameters_in.insert_name_value("w2", 1.2);

    if (os_) *os_ << "parameters_in: " << endl << parameters_in << endl;
    
    QuantitativeTrait_SingleLocusFitness qt("my_qt");
    qt.configure(parameters_in, registry);

    Parameters parameters_out = qt.parameters();
    if (os_) *os_ << "parameters_in: " << endl << parameters_in << endl;

    unit_assert(parameters_in == parameters_out);

    set<string> ids_written;
    if (os_)
    {
        *os_ << "write_configuration():\n\n";
        qt.write_configuration(*os_, ids_written);
    }
}


class QuantitativeTrait_TestingComposite : public QuantitativeTrait
{
    public:

    QuantitativeTrait_TestingComposite(const string& id) : QuantitativeTrait(id) {}

    virtual void calculate_trait_values(const PopulationData& population_data) const
    {
        population_indices.push_back(population_data.population_index);
    }

    mutable vector<size_t> population_indices;
};


typedef shared_ptr<QuantitativeTrait_TestingComposite> QuantitativeTrait_TestingCompositePtr;

 
void test_QuantitativeTrait_PopulationComposite()
{
    if (os_) *os_ << "test_QuantitativeTrait_PopulationComposite()\n";

    string qtid_testing_0 = "qt_testing_0";
    string qtid_testing_1 = "qt_testing_1";
    string qtid_testing_2 = "qt_testing_2";

    QuantitativeTraitPtr qt_testing_0(new QuantitativeTrait_TestingComposite(qtid_testing_0));
    QuantitativeTraitPtr qt_testing_1(new QuantitativeTrait_TestingComposite(qtid_testing_1));
    QuantitativeTraitPtr qt_testing_2(new QuantitativeTrait_TestingComposite(qtid_testing_2));

    Parameters parameters;
    ostringstream qts_string;
    qts_string << qtid_testing_0 << " " << qtid_testing_1 << " " << qtid_testing_2 << " ";
    parameters.insert_name_value("quantitative_traits", qts_string.str());

    Configurable::Registry registry;
    registry[qtid_testing_0] = qt_testing_0;
    registry[qtid_testing_1] = qt_testing_1;
    registry[qtid_testing_2] = qt_testing_2;

    QuantitativeTrait_PopulationComposite qt("qt_population_composite");
    qt.configure(parameters, registry);

    Parameters parameters_out = qt.parameters();

    set<string> ids_written;
    if (os_) qt.write_configuration(*os_, ids_written);

    unit_assert(parameters == parameters_out);

    PopulationData population_data;

    population_data.generation_index = 0;
    population_data.population_index = 0;
    qt.calculate_trait_values(population_data);

    population_data.population_index = 1;
    qt.calculate_trait_values(population_data);
    population_data.trait_values->clear();
    qt.calculate_trait_values(population_data);

    population_data.population_index = 2;
    qt.calculate_trait_values(population_data);
    population_data.trait_values->clear();
    qt.calculate_trait_values(population_data);
    population_data.trait_values->clear();
    qt.calculate_trait_values(population_data);

    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_0)
            ->population_indices.size() == 1);
    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_1)
            ->population_indices.size() == 2);
    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_2)
            ->population_indices.size() == 3);
}


void test_QuantitativeTrait_GenerationComposite()
{
    if (os_) *os_ << "test_QuantitativeTrait_GenerationComposite()\n";

    string qtid_testing_0 = "qt_testing_0";
    string qtid_testing_1 = "qt_testing_1";
    string qtid_testing_2 = "qt_testing_2";

    QuantitativeTraitPtr qt_testing_0(new QuantitativeTrait_TestingComposite(qtid_testing_0));
    QuantitativeTraitPtr qt_testing_1(new QuantitativeTrait_TestingComposite(qtid_testing_1));
    QuantitativeTraitPtr qt_testing_2(new QuantitativeTrait_TestingComposite(qtid_testing_2));

    Parameters parameters;
    parameters.insert_name_value("generation:quantitative_trait", "0 " + qtid_testing_0);
    parameters.insert_name_value("generation:quantitative_trait", "5 " + qtid_testing_1);
    parameters.insert_name_value("generation:quantitative_trait", "23 " + qtid_testing_2);

    Configurable::Registry registry;
    registry[qtid_testing_0] = qt_testing_0;
    registry[qtid_testing_1] = qt_testing_1;
    registry[qtid_testing_2] = qt_testing_2;

    QuantitativeTrait_GenerationComposite qt("qt_population_composite");
    qt.configure(parameters, registry);

    Parameters parameters_out = qt.parameters();

    set<string> ids_written;
    if (os_) qt.write_configuration(*os_, ids_written);

    unit_assert(parameters == parameters_out);

    PopulationData population_data;

    for (size_t generation_index=0; generation_index<30; ++generation_index)
    {
        population_data.generation_index = generation_index;
        population_data.trait_values->clear();
        qt.calculate_trait_values(population_data);
    }

    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_0)
            ->population_indices.size() == 5);
    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_1)
            ->population_indices.size() == 18);
    unit_assert(dynamic_pointer_cast<QuantitativeTrait_TestingComposite>(qt_testing_2)
            ->population_indices.size() == 7);
}


void test_QuantitativeTrait_IndependentLoci()
{
    if (os_) *os_ << "test_QuantitativeTrait_IndependentLoci()\n";

    LocusPtr locus_1(new Locus("id_locus_1", 0, 1000));
    LocusPtr locus_2(new Locus("id_locus_2", 0, 2000));
    LocusPtr locus_3(new Locus("id_locus_3", 0, 3000));

    QTLEffects qtl_effects;
    qtl_effects.push_back(QTLEffect(*locus_1, 0, 100, 200));
    qtl_effects.push_back(QTLEffect(*locus_2, 0, 10, 20));
    qtl_effects.push_back(QTLEffect(*locus_3, 0, 1, 2));

    QuantitativeTrait_IndependentLoci qt(string("id_dummy"), qtl_effects);

    // individual genotypes: locus_1 locus_2 locus_3  trait_value
    //      0 0 0  000
    //      1 1 1  111
    //      2 2 2  222
    //      0 1 2  012
    //      2 1 0  210

    size_t n = 5;
    char genotypes_raw_1[5] = {0, 1, genotype_make_pair(1,1), 0, genotype_make_pair(1,1)};
    char genotypes_raw_2[5] = {0, 1, genotype_make_pair(1,1), 1, 1};
    char genotypes_raw_3[5] = {0, 1, genotype_make_pair(1,1), genotype_make_pair(1,1), 0};

    PopulationData population_data;
    population_data.generation_index = 0;
    population_data.population_index = 0;
    population_data.population_size = n;

    GenotypeMapPtr genotypes = population_data.genotypes;
    (*genotypes)[*locus_1] = GenotypeDataPtr(new GenotypeData(genotypes_raw_1, genotypes_raw_1 + n));
    (*genotypes)[*locus_2] = GenotypeDataPtr(new GenotypeData(genotypes_raw_2, genotypes_raw_2 + n));
    (*genotypes)[*locus_3] = GenotypeDataPtr(new GenotypeData(genotypes_raw_3, genotypes_raw_3 + n));

    qt.calculate_trait_values(population_data);
    DataVectorPtr trait_values = population_data.trait_values->at("id_dummy");
    unit_assert(trait_values.get());

    if (os_)
        *os_ << "trait_values: " << *trait_values << endl;

    unit_assert(trait_values->size() == n);
    unit_assert((*trait_values)[0] == 0);
    unit_assert((*trait_values)[1] == 111);
    unit_assert((*trait_values)[2] == 222);
    unit_assert((*trait_values)[3] == 12);
    unit_assert((*trait_values)[4] == 210);

    // re-configure to nothing

    Parameters parameters;

    Configurable::Registry registry;
    registry["id_locus_1"] = locus_1;
    registry["id_locus_2"] = locus_2;
    registry["id_locus_3"] = locus_3;

    qt.configure(parameters, registry);

    PopulationData population_data_empty; // no genotypes
    population_data_empty.population_size = n;

    qt.calculate_trait_values(population_data_empty);
    trait_values = population_data_empty.trait_values->at("id_dummy");
    unit_assert(trait_values.get());
    unit_assert(trait_values->size() == n);

    if (os_) *os_ << "trait_values: " << *trait_values << endl;

    for (DataVector::const_iterator it=trait_values->begin(); it!=trait_values->end(); ++it)
        unit_assert(*it == 0);

    // configure: environmental variance only

    parameters.insert_name_value("environmental_variance", 1);
    qt.configure(parameters, registry);

    SimulatorConfig simconfig;
    qt.initialize(simconfig); // secondary initialization for environmental variance distribution

    population_data_empty.trait_values->clear(); // reset trait values
    qt.calculate_trait_values(population_data_empty);
    trait_values = population_data_empty.trait_values->at("id_dummy");

    unit_assert(trait_values.get());
    unit_assert(trait_values->size() == n);

    if (os_) *os_ << "trait_values: " << *trait_values << endl;
    for (DataVector::const_iterator it=trait_values->begin(); it!=trait_values->end(); ++it)
        unit_assert(*it != 0);

    // configure: original qtls

    parameters.clear();
    parameters.insert_name_value("qtl", "id_locus_1 0 100 200");
    parameters.insert_name_value("qtl", "id_locus_2 0 10 20");
    parameters.insert_name_value("qtl", "id_locus_3 0 1 2");
    parameters.insert_name_value("environmental_variance", 0);

    qt.configure(parameters, registry);
    qt.initialize(simconfig); // secondary initialization for environmental variance distribution

    if (os_)
    {
        set<string> written;
        qt.write_configuration(*os_, written);
    }

    Parameters parameters_out = qt.parameters();
    unit_assert(parameters == parameters_out);

    population_data.trait_values->clear(); // reset trait values
    qt.calculate_trait_values(population_data);
    trait_values = population_data.trait_values->at("id_dummy");

    if (os_) *os_ << "trait_values: " << *trait_values << endl;

    unit_assert(trait_values.get());
    unit_assert(trait_values->size() == n);
    unit_assert((*trait_values)[0] == 0);
    unit_assert((*trait_values)[1] == 111);
    unit_assert((*trait_values)[2] == 222);
    unit_assert((*trait_values)[3] == 12);
    unit_assert((*trait_values)[4] == 210);
}


void test_QTLEffectGenerator()
{
    if (os_) *os_ << "test_QTLEffectGenerator()\n";

    Random::seed(1);

    LocusListPtr locus_list(new LocusList("locus_list"));
    locus_list->push_back(Locus("locus_list_1", 1, 2345));
    locus_list->push_back(Locus("locus_list_2", 2, 3456));
    locus_list->push_back(Locus("locus_list_3", 3, 4567));

    Random::DistributionPtr effect_size_distribution = 
        Random::create_constant_distribution("effect_size_distribution", 1);

    Random::DistributionPtr dominance_distribution = 
        Random::create_constant_distribution("dominance_distribution", 0);

    QTLEffectGeneratorPtr qtl_effect_generator(new QTLEffectGenerator("qtl_effect_generator",
                                                                      locus_list,
                                                                      effect_size_distribution,
                                                                      dominance_distribution));

    QTLEffects qtl_effects;

    qtl_effect_generator->generate_qtl_effects(qtl_effects);

    unit_assert(qtl_effects.size() == locus_list->size());

    // plug into QT

    Random::seed(1);

    QuantitativeTrait_IndependentLoci qt(string("id_dummy"));

    LocusPtr locus_1(new Locus("id_locus_1", 0, 1000));
    LocusPtr locus_2(new Locus("id_locus_2", 0, 2000));
    LocusPtr locus_3(new Locus("id_locus_3", 0, 3000));

    Configurable::Registry registry;
    registry["qtl_effect_generator"] = qtl_effect_generator;
    registry["id_locus_1"] = locus_1;
    registry["id_locus_2"] = locus_2;
    registry["id_locus_3"] = locus_3;

    Parameters parameters;
    parameters.insert_name_value("qtl", "id_locus_1 0 100 200");
    parameters.insert_name_value("qtl", "id_locus_2 0 10 20");
    parameters.insert_name_value("qtl", "id_locus_3 0 1 2");
    parameters.insert_name_value("qtl_effect_generator", "qtl_effect_generator");
    parameters.insert_name_value("environmental_variance", 0);

    qt.configure(parameters, registry);

    SimulatorConfig simconfig;
    qt.initialize(simconfig);

    const QTLEffects& qtl_effects_out = qt.qtl_effects();

    if (os_)
    {
        *os_ << "qtl_effects_out:\n";
        for (QTLEffects::const_iterator it=qtl_effects_out.begin(); it!=qtl_effects_out.end(); ++it)
            *os_ << it->locus << " " << it->configuration() << endl;
        *os_ << endl;
    }

    unit_assert(qtl_effects_out.size() == locus_list->size() + 3);
    const Loci& loci = qt.loci();
    unit_assert(loci.size() == locus_list->size() + 3);

    // check configuration

    Parameters parameters_out = qt.parameters();

    if (os_)
    {
        set<string> written;
        qt.write_configuration(*os_, written);

        *os_ << "parameters:\n" << parameters << endl
             << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_Configurable_QTLEffectGenerator()
{
    if (os_) *os_ << "test_Configurable_QTLEffectGenerator()\n" << flush;

    LocusListPtr locus_list = LocusListPtr(new LocusList("locus_list"));

    Random::DistributionPtr effect_size_distribution = 
        Random::create_constant_distribution("id_effect_size_distribution", 1);

    Random::DistributionPtr dominance_distribution = 
        Random::create_constant_distribution("id_dominance_distribution", 0);

    Configurable::Registry registry;
    registry[locus_list->object_id()] = locus_list;
    registry[effect_size_distribution->object_id()] = effect_size_distribution;
    registry[dominance_distribution->object_id()] = dominance_distribution;

    Parameters parameters;
    parameters.insert_name_value("locus_list", locus_list->object_id());
    parameters.insert_name_value("effect_size_distribution", effect_size_distribution->object_id());
    parameters.insert_name_value("dominance_distribution", dominance_distribution->object_id());

    QTLEffectGeneratorPtr qtl_effect_generator(new QTLEffectGenerator("qtl_effect_generator"));

    // configure and check
    
    qtl_effect_generator->configure(parameters, registry);

    Parameters parameters_out = qtl_effect_generator->parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
        *os_ << "configuration:\n";
        set<string> written;
        qtl_effect_generator->write_configuration(*os_, written);
    }

    unit_assert(parameters == parameters_out);
}


void test_QuantitativeTrait_Expression()
{
    if (os_) *os_ << "test_QuantitativeTrait_Expression()\n";

    const size_t sample_size = 10;

    DataVectorPtr index(new DataVector(sample_size));
    DataVectorPtr even(new DataVector(sample_size));

    for (size_t i=0; i<sample_size; ++i)
    {
        (*index)[i] = i;
        (*even)[i] = (i%2==0);
    }

    if (os_)
    {
        *os_ << "index: " << *index << endl;
        *os_ << "even: " << *even << endl;
    }

    PopulationData popdata;
    popdata.population_size = sample_size;
    TraitValueMap& trait_values = *popdata.trait_values;
    trait_values["index"] = index;
    trait_values["even"] = even;

    Parameters parameters_in;
    parameters_in.insert_name_value("variable:quantitative_trait", "x index");
    parameters_in.insert_name_value("variable:quantitative_trait", "y even");
    parameters_in.insert_name_value("expression", "x * y");

    Configurable::Registry registry;
    QuantitativeTrait_Expression qt("qt_expression");
    qt.configure(parameters_in, registry);

    Parameters parameters_out = qt.parameters();

    if (os_)
    {
        *os_ << "parameters_in: " << endl << parameters_in << endl;
        *os_ << "parameters_out: " << endl << parameters_out << endl;
    }

    unit_assert(parameters_in == parameters_out);

    // test 1: x * y

    qt.calculate_trait_values(popdata);
    unit_assert(trait_values.size() == 3 && trait_values.count("qt_expression"));
    const DataVector& values = *trait_values.get("qt_expression");
    if (os_)
    {
        *os_ << "x: " << *index << endl
             << "y: " << *even << endl
             << "x*y: " << values << endl;
    }
    for (size_t i=0; i<10; i+=2) unit_assert(values.at(i) == i);
    for (size_t i=1; i<10; i+=2) unit_assert(values.at(i) == 0);

    // test 2: 2^y * x

    parameters_in.erase("expression");
    trait_values.erase("qt_expression");

    parameters_in.insert_name_value("expression", "x * 2^y");
    qt.configure(parameters_in, registry);
    qt.calculate_trait_values(popdata);

    const DataVector& values2 = *trait_values.get("qt_expression");
    if (os_) *os_ << "x * 2^y: " << values2 << endl;
    for (size_t i=0; i<10; i+=2) unit_assert(values2.at(i) == 2*i);
    for (size_t i=1; i<10; i+=2) unit_assert(values2.at(i) == i);

    // test 3: (x<6)*y

    parameters_in.erase("expression");
    trait_values.erase("qt_expression");

    parameters_in.insert_name_value("expression", "(x<6)*y");
    qt.configure(parameters_in, registry);
    qt.calculate_trait_values(popdata);

    const DataVector& values3 = *trait_values.get("qt_expression");
    if (os_) *os_ << "(x<6)*y: " << values3 << endl;
    for (size_t i=0; i<6; ++i) unit_assert(values3.at(i) == (i%2==0));
    for (size_t i=6; i<10; ++i) unit_assert(values3.at(i) == 0);
}


void test()
{
    test_Configurable_QuantitativeTrait_SingleLocusFitness();
    test_QuantitativeTrait_PopulationComposite();
    test_QuantitativeTrait_GenerationComposite();
    test_QuantitativeTrait_IndependentLoci();
    test_QTLEffectGenerator();
    test_Configurable_QTLEffectGenerator();
    test_QuantitativeTrait_Expression();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        test();
        return 0;
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


