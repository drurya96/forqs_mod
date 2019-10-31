//
// ReporterImplementationTest.cpp
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


#include "ReporterImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <numeric>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_Reporter_Timer()
{
    const size_t update_count = 3;

    Reporter_Timer reporter("reporter_timer");
    for (size_t i=0; i<update_count-1; i++)
    {
        time_t begin = time(0);
        while (difftime(time(0), begin) < 1); // wait 1 second
        reporter.update(0, PopulationPtrs(), PopulationDataPtrs(), false);
    }
    reporter.update(0, PopulationPtrs(), PopulationDataPtrs(), true);

    const vector<double>& times = reporter.times();

    // verify number of updates

    unit_assert(times.size() == update_count);
    
    // verify that inner updates are close to 1 second intervals

    vector<double> diffs;
    adjacent_difference(times.begin(), times.end(), back_inserter(diffs));

    //copy(times.begin(), times.end(), ostream_iterator<double>(cout, " ")); cout << endl;
    //copy(diffs.begin(), diffs.end(), ostream_iterator<double>(cout, " ")); cout << endl;

    const double epsilon = .11;
    for (size_t i=1; i<diffs.size()-1; ++i)
        unit_assert_equal(diffs[i], 1.0, epsilon);

    //cout << "mean time: " << reporter.mean_generation_time() << endl;
    unit_assert_equal(reporter.mean_generation_time(), 1.0, epsilon);
}


void test_Configurable_Reporter_Population()
{
    if (os_) *os_ << "test_Configurable_Reporter_Population()\n";

    Parameters parameters_in;
    parameters_in.insert_name_value("update_step", "10");

    if (os_) *os_ << "parameters_in:\n" << parameters_in << endl;

    Reporter_Population reporter("my_reporter");

    Configurable::Registry registry;
    reporter.configure(parameters_in, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_) *os_ << "parameters_out:\n" << parameters_out << endl;

    unit_assert(parameters_in == parameters_out);
}


namespace {
class QuantitativeTrait_Testing : public QuantitativeTrait
{
    public:

    QuantitativeTrait_Testing()
    :   QuantitativeTrait("qt_testing")
    {
        loci_.insert(Locus("qt_1", 3, 456789));
        loci_.insert(Locus("qt_1", 4, 567890));
    }

    const Loci& loci() const {return loci_;}

    virtual void calculate_trait_values(const PopulationData& population_data) const {}
};
} // namespace


void test_Configurable_Reporter_AlleleFrequencies()
{
    if (os_) *os_ << "test_Configurable_Reporter_AlleleFrequencies()\n";

    LocusPtr locus(new Locus("my_locus", 0, 123456));

    LocusListPtr locus_list(new LocusList("my_locus_list"));
    locus_list->push_back(Locus("locus_list_1", 1, 234567));
    locus_list->push_back(Locus("locus_list_2", 2, 345678));

    QuantitativeTraitPtr qt = QuantitativeTraitPtr(new QuantitativeTrait_Testing);

    Configurable::Registry registry;
    registry[locus->object_id()] = locus;
    registry[locus_list->object_id()] = locus_list;
    registry[qt->object_id()] = qt;

    Parameters parameters_in;
    parameters_in.insert_name_value("locus", locus->object_id());
    parameters_in.insert_name_value("locus_list", locus_list->object_id());
    parameters_in.insert_name_value("quantitative_trait", qt->object_id());

    Reporter_AlleleFrequencies reporter("dummy_id");
    reporter.configure(parameters_in, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "locus_list size: " << locus_list->size() << endl;
        *os_ << "qt size: " << qt->loci().size() << "\n\n";

        *os_ << "parameters_in:\n" << parameters_in << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "write_configuration():\n\n";
        set<string> ids_written;
        reporter.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters_in == parameters_out);
}


void test_Configurable_Reporter_LD()
{
    if (os_) *os_ << "test_Configurable_Reporter_LD()\n";

    string id_locus1 = "id_locus1";
    string id_locus2 = "id_locus2";

    LocusPtr locus1(new Locus(id_locus1, 0, 123456));
    LocusPtr locus2(new Locus(id_locus2, 1, 789012));

    Configurable::Registry registry;
    registry[id_locus1] = locus1;
    registry[id_locus2] = locus2;

    Parameters parameters_in;
    parameters_in.insert_name_value("locus_1", id_locus1);
    parameters_in.insert_name_value("locus_2", id_locus2);

    Reporter_LD reporter("dummy_id");
    reporter.configure(parameters_in, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_) 
    {
        *os_ << "parameters_in:\n" << parameters_in << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        reporter.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters_in == parameters_out);
}


void test_Configurable_Reporter_TraitValues()
{
    if (os_) *os_ << "test_Configurable_Reporter_TraitValues()\n";

    Parameters parameters;
    parameters.insert_name_value("quantitative_traits", "qt1 qt2 qt3 ");
    parameters.insert_name_value("ignore_zero_values", "qt2 ");
    parameters.insert_name_value("write_full", 1);
    parameters.insert_name_value("filetag", "goober");

    Configurable::Registry registry;

    Reporter_TraitValues reporter("trait_values");

    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        reporter.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_Configurable_Reporter_HaplotypeDiversity()
{
    if (os_) *os_ << "test_Configurable_Reporter_HaplotypeDiversity()\n";

    Parameters parameters;
    parameters.insert_name_value("chromosome", "0 1000000 10000");
    parameters.insert_name_value("chromosome", "3 2000000 1000");

    Configurable::Registry registry;

    Reporter_HaplotypeDiversity reporter("id_dummy");

    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        reporter.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_Configurable_Reporter_HaplotypeFrequencies()
{
    if (os_) *os_ << "test_Configurable_Reporter_HaplotypeFrequencies()\n";

    HaplotypeGroupingPtr haplotype_grouping(new HaplotypeGrouping_IDRange("id_grouping"));

    Configurable::Registry registry;
    registry["id_grouping"] = haplotype_grouping;

    Parameters parameters;
    parameters.insert_name_value("haplotype_grouping", "id_grouping");
    parameters.insert_name_value("chromosome_step", 100000);
    parameters.insert_name_value("update_step", 100);
    
    Reporter_HaplotypeFrequencies reporter("id_dummy");
    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        reporter.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_HaplotypeGrouping_IDRange()
{
    if (os_) *os_ << "test_HaplotypeGrouping_IDRange()\n";

    Parameters parameters;
    parameters.insert_name_value("start:count", "0 100");
    parameters.insert_name_value("start:count", "1000 100");
    parameters.insert_name_value("start:count", "2000 100");

    HaplotypeGrouping_IDRange grouping("id_dummy");

    Configurable::Registry registry;
    grouping.configure(parameters, registry);

    Parameters parameters_out = grouping.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    unit_assert(grouping.group_count() == 3);
    unit_assert(grouping.group(0) == 0);
    unit_assert(grouping.group(1) == 0);
    unit_assert(grouping.group(10) == 0);
    unit_assert(grouping.group(90) == 0);
    unit_assert(grouping.group(99) == 0);
    unit_assert(grouping.group(1000) == 1);
    unit_assert(grouping.group(1001) == 1);
    unit_assert(grouping.group(1099) == 1);
    unit_assert(grouping.group(2000) == 2);
    unit_assert(grouping.group(2002) == 2);
    unit_assert(grouping.group(2099) == 2);
}


class PCG_Test : public PopulationConfigGenerator
{
    public:

    PCG_Test() : PopulationConfigGenerator("dummy") {}

    virtual size_t generation_count() const {return 1;}

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const
    {
        Population::Configs configs(1);
        configs[0].id_offset = 3;
        configs[0].population_size = 1000;
        return configs;
    }
};


void test_HaplotypeGrouping_Uniform()
{
    if (os_) *os_ << "test_HaplotypeGrouping_Uniform()\n";

    Parameters parameters;
    parameters.insert_name_value("ids_per_group", 2);

    HaplotypeGrouping_Uniform grouping("id_dummy");

    Configurable::Registry registry;
    grouping.configure(parameters, registry);

    Parameters parameters_out = grouping.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    // secondary initialization
    boost::shared_ptr<PCG_Test> pcg(new PCG_Test);
    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;
    unit_assert(pcg->population_configs(0,PopulationDataPtrs())[0].population_size == 1000);
    const unsigned int population_size = pcg->population_configs(0,PopulationDataPtrs())[0].population_size;
    const unsigned int id_offset = pcg->population_configs(0,PopulationDataPtrs())[0].id_offset;
    grouping.initialize(simconfig);

    if (os_)
        *os_ << "group count: " << grouping.group_count() << endl;

    unit_assert(grouping.group_count() == 1000); // N=1000 * 2 ids/indiv == 2000 ids total

    for (unsigned int i=0; i<population_size; ++i)
        unit_assert(grouping.group(i) == (i-id_offset)/2);
}


void test_Configurable_Reporter_Regions()
{
    if (os_) *os_ << "test_Reporter_Regions()\n";

    LocusPtr locus1(new Locus("locus1", 0, 123456));
    LocusPtr locus2(new Locus("locus2", 0, 789012));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;

    Parameters parameters;
    parameters.insert_name_value("locus:length", "locus1 1000");
    parameters.insert_name_value("locus:length", "locus2 2000");
    parameters.insert_name_value("ms_mapping_begin", 1234);
    parameters.insert_name_value("ms_mapping_end", 5678);

    Reporter_Regions reporter("id_dummy");
    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        reporter.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_Configurable_Reporter_DeterministicTrajectories()
{
    if (os_) *os_ << "test_Reporter_DeterministicTrajectories()\n";

    Parameters parameters;
    parameters.insert_name_value("initial_allele_frequency", .1);
    parameters.insert_name_value("w0", 1);
    parameters.insert_name_value("w1", 1.1);
    parameters.insert_name_value("w2", 1.2);

    Configurable::Registry registry;

    Reporter_DeterministicTrajectories reporter("id_dummy");
    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        reporter.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);
}


void test_Configurable_Reporter_Variants()
{
    if (os_) *os_ << "test_Configurable_Reporter_Variants()\n";

    LocusPtr locus1(new Locus("locus1", 1, 20000));
    LocusPtr locus2(new Locus("locus2", 2, 40000));
    LocusPtr locus3(new Locus("locus3", 3, 60000));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;
    registry["locus3"] = locus3;

    LocusListPtr locus_list(new LocusList("locus_list"));
    Parameters locus_list_parameters;
    locus_list_parameters.insert_name_value("loci", "locus2 locus3");
    locus_list->configure(locus_list_parameters, registry);
    registry["locus_list"] = locus_list;

    QuantitativeTraitPtr qt = QuantitativeTraitPtr(new QuantitativeTrait_Testing);
    registry[qt->object_id()] = qt;

    Parameters parameters;
    parameters.insert_name_value("loci", "locus1 locus_list qt_testing ");
    parameters.insert_name_value("update_step", 3);

    Reporter_Variants reporter("id_dummy");
    reporter.configure(parameters, registry);

    Parameters parameters_out = reporter.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        reporter.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);
}


void test()
{
    test_Reporter_Timer();
    test_Configurable_Reporter_Population();
    test_Configurable_Reporter_AlleleFrequencies();
    test_Configurable_Reporter_LD();
    test_Configurable_Reporter_TraitValues();
    test_Configurable_Reporter_HaplotypeDiversity();
    test_Configurable_Reporter_HaplotypeFrequencies();
    test_HaplotypeGrouping_IDRange();
    test_HaplotypeGrouping_Uniform();
    test_Configurable_Reporter_Regions();
    test_Configurable_Reporter_DeterministicTrajectories();
    test_Configurable_Reporter_Variants();
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


