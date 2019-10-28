//
// VariantIndicatorImplementationTest.cpp
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


#include "VariantIndicatorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


class PCG_TestHW : public PopulationConfigGenerator
{
    public:

    PCG_TestHW() : PopulationConfigGenerator("dummy") {}

    virtual size_t generation_count() const {return 1;}

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const
    {
        Population::Configs configs(3);
        for (size_t i=0; i<3; ++i)
        {
            configs[i].id_offset = i*10000;
            configs[i].population_size = (i+1)*1000;
        }
        return configs;
    }
};


void test_VariantIndicator_SingleLocusHardyWeinberg()
{
    if (os_) *os_ << "test_VariantIndicator_SingleLocusHardyWeinberg()\n";

    Locus locus("id_dummy");
    const double allele_frequency = .5;
    VariantIndicator_SingleLocusHardyWeinberg vi("dummy_id", locus, allele_frequency);

    shared_ptr<PCG_TestHW> pcg(new PCG_TestHW);
    Population::Configs popconfigs = pcg->population_configs(0, PopulationDataPtrs());

    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    vi.initialize(simconfig);

    for (Population::Configs::const_iterator popconfig=popconfigs.begin(); popconfig!=popconfigs.end(); ++popconfig)
    {
        size_t start = popconfig->id_offset;
        size_t count = 2 * popconfig->population_size;

        if (os_) *os_ << "ids: " << start << " " << count << endl;

        for (size_t id = start; id < start + count/4; ++id)
            unit_assert(vi(id,locus) == 1);

        for (size_t id = start + count/4; id < start + count*3/4; ++id)
            unit_assert(vi(id,locus) == ((id%2==0) ? 1 : 0));

        for (size_t id = start+count*3/4; id < start+count; ++id)
            unit_assert(vi(id,locus) == 0);
    }
}


void test_Configurable_VariantIndicator_SingleLocusHardyWeinberg()
{
    if (os_) *os_ << "test_Configurable_VariantIndicator_SingleLocusHardyWeinberg()\n";

    string locus_id = "my_locus";
    LocusPtr locus(new Locus(locus_id, 0, 123456));

    Configurable::Registry registry;
    registry[locus_id] = locus;

    Parameters parameters_in;
    parameters_in.insert_name_value("locus", locus_id);
    parameters_in.insert_name_value("allele_frequency", .25);

    if (os_) *os_ << "parameters_in:\n" << parameters_in << endl;

    VariantIndicatorPtr snpi(new VariantIndicator_SingleLocusHardyWeinberg("dummy_id"));
    snpi->configure(parameters_in, registry);

    Parameters parameters_out = snpi->parameters();

    if (os_) *os_ << "parameters_out:\n" << parameters_out << endl;

    unit_assert(parameters_in == parameters_out);
    
    set<string> ids_written;
    if (os_) 
    {
        *os_ << "write_configuration():\n\n";
        snpi->write_configuration(*os_, ids_written);
    }
}


void test_VariantIndicator_TwoLocusLD()
{
    if (os_) *os_ << "test_VariantIndicator_TwoLocusLD()\n";

    const size_t N = 100;
    const double p00 = .4;
    const double p01 = .1;
    const double p10 = .1;
    const double p11 = .4;
    const double D = p00*p11 - p01*p10;
    const size_t id_offset_step = 1000;

    string id_locus_1 = "id_locus_1";
    string id_locus_2 = "id_locus_2";
    LocusPtr locus1(new Locus(id_locus_1, 0, 100000));
    LocusPtr locus2(new Locus(id_locus_2, 0, 200000));

    Configurable::Registry registry;
    registry[id_locus_1] = locus1;
    registry[id_locus_2] = locus2;

    Parameters parameters;
    parameters.insert_name_value("population_size", N);
    parameters.insert_name_value("locus_1", id_locus_1);
    parameters.insert_name_value("allele_frequency_1", p10 + p11);
    parameters.insert_name_value("locus_2", id_locus_2);
    parameters.insert_name_value("allele_frequency_2", p01 + p11);
    parameters.insert_name_value("D", D);
    parameters.insert_name_value("id_offset_step", id_offset_step);

    VariantIndicator_TwoLocusLD s("id_dummy");
    s.configure(parameters, registry);

    Parameters parameters_out = s.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> written;
        s.write_configuration(*os_, written);
    }

    unit_assert(parameters == parameters_out);

    vector<unsigned int> snp_values_1;
    for (size_t i=0; i<N; ++i) snp_values_1.push_back(0);
    for (size_t i=0; i<N; ++i) snp_values_1.push_back(1);
    
    vector<unsigned int> snp_values_2;
    for (size_t i=0; i<2*N*p00; ++i) snp_values_2.push_back(0);
    for (size_t i=0; i<2*N*p01; ++i) snp_values_2.push_back(1);
    for (size_t i=0; i<2*N*p10; ++i) snp_values_2.push_back(0);
    for (size_t i=0; i<2*N*p11; ++i) snp_values_2.push_back(1);
    unit_assert(snp_values_2.size() == 2*N);

    for (size_t i=0; i<2*N; ++i)
    {
        unit_assert(s(i, *locus1) == snp_values_1[i]);
        unit_assert(s(i, *locus2) == snp_values_2[i]);
    }
}


void test_VariantIndicator_IDRange()
{
    if (os_) *os_ << "test_VariantIndicator_IDRange()\n";
    
    LocusPtr locus1(new Locus("locus1", 1, 10000));
    LocusPtr locus2(new Locus("locus2", 2, 20000));
    LocusPtr locus3(new Locus("locus3", 3, 30000));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;
    registry["locus3"] = locus3;

    Parameters parameters;
    parameters.insert_name_value("locus:start:count:step:value", "locus1 0 1000 1 10");
    parameters.insert_name_value("locus:start:count:step:value", "locus2 1000 1000 2 11");
    parameters.insert_name_value("locus:start:count:step:value", "locus3 2000 1000 3 12");

    VariantIndicator_IDRange vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    unit_assert(vi(0,*locus1) == 10);
    unit_assert(vi(1,*locus1) == 10);
    unit_assert(vi(2,*locus1) == 10);
    unit_assert(vi(999,*locus1) == 10);
    unit_assert(vi(1000,*locus1) == 0);
    unit_assert(vi(2000,*locus1) == 0);
    unit_assert(vi(3000,*locus1) == 0);

    unit_assert(vi(1000,*locus2) == 11);
    unit_assert(vi(1001,*locus2) == 0);
    unit_assert(vi(1002,*locus2) == 11);
    unit_assert(vi(1003,*locus2) == 0);
    unit_assert(vi(1998,*locus2) == 11);
    unit_assert(vi(1999,*locus2) == 0);
    unit_assert(vi(0,*locus2) == 0);
    unit_assert(vi(2000,*locus2) == 0);
    unit_assert(vi(3000,*locus2) == 0);

    unit_assert(vi(2000,*locus3) == 12);
    unit_assert(vi(2001,*locus3) == 0);
    unit_assert(vi(2002,*locus3) == 0);
    unit_assert(vi(2003,*locus3) == 12);
    unit_assert(vi(2997,*locus3) == 0);
    unit_assert(vi(2998,*locus3) == 0);
    unit_assert(vi(2999,*locus3) == 12);
    unit_assert(vi(0,*locus3) == 0);
    unit_assert(vi(1000,*locus3) == 0);
    unit_assert(vi(3000,*locus3) == 0);

    Locus locus4("locus4", 0, 100000);
    unit_assert(vi(0,locus4) == 0);
    unit_assert(vi(1000,locus4) == 0);
    unit_assert(vi(2000,locus4) == 0);
    unit_assert(vi(3000,locus4) == 0);
}


void test_VariantIndicator_IDSet()
{
    if (os_) *os_ << "test_VariantIndicator_IDSet()\n";
    
    LocusPtr locus1(new Locus("locus1", 1, 10000));
    LocusPtr locus2(new Locus("locus2", 2, 20000));
    LocusPtr locus3(new Locus("locus3", 3, 30000));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;
    registry["locus3"] = locus3;

    Parameters parameters;
    parameters.insert_name_value("locus:value:ids", "locus1 11 1001 1002 1003 ");
    parameters.insert_name_value("locus:value:ids", "locus2 22 2001 2002 2003 2004 ");
    parameters.insert_name_value("locus:value:ids", "locus3 33 3001 3002 3003 3004 3005 ");

    VariantIndicator_IDSet vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    unit_assert(vi(0,*locus1) == 0);
    unit_assert(vi(1,*locus1) == 0);
    unit_assert(vi(1000,*locus1) == 0);
    unit_assert(vi(1001,*locus1) == 11);
    unit_assert(vi(1002,*locus1) == 11);
    unit_assert(vi(1003,*locus1) == 11);
    unit_assert(vi(1004,*locus1) == 0);

    unit_assert(vi(0,*locus2) == 0);
    unit_assert(vi(1,*locus2) == 0);
    unit_assert(vi(2000,*locus2) == 0);
    unit_assert(vi(2001,*locus2) == 22);
    unit_assert(vi(2002,*locus2) == 22);
    unit_assert(vi(2003,*locus2) == 22);
    unit_assert(vi(2004,*locus2) == 22);
    unit_assert(vi(2005,*locus2) == 0);

    unit_assert(vi(0,*locus3) == 0);
    unit_assert(vi(1,*locus3) == 0);
    unit_assert(vi(3000,*locus3) == 0);
    unit_assert(vi(3001,*locus3) == 33);
    unit_assert(vi(3002,*locus3) == 33);
    unit_assert(vi(3003,*locus3) == 33);
    unit_assert(vi(3004,*locus3) == 33);
    unit_assert(vi(3005,*locus3) == 33);
    unit_assert(vi(3006,*locus3) == 0);

    Locus locus4("locus4", 4, 40000);
    unit_assert(vi(0,locus4) == 0);
    unit_assert(vi(1,locus4) == 0);
    unit_assert(vi(1000,locus4) == 0);
    unit_assert(vi(1001,locus4) == 0);
    unit_assert(vi(2000,locus4) == 0);
    unit_assert(vi(2001,locus4) == 0);
    unit_assert(vi(3000,locus4) == 0);
    unit_assert(vi(3001,locus4) == 0);
}


void test_VariantIndicator_Random()
{
    if (os_) *os_ << "test_VariantIndicator_Random()\n\n";
    
    Random::seed(2);

    LocusPtr locus1(new Locus("locus1", 1, 10000));
    LocusPtr locus2(new Locus("locus2", 2, 20000));
    LocusPtr locus3(new Locus("locus3", 3, 30000));
    LocusPtr locus4(new Locus("locus4", 4, 40000));

    Configurable::Registry registry;
    registry[locus1->object_id()] = locus1;
    registry[locus2->object_id()] = locus2;
    registry[locus3->object_id()] = locus3;
    registry[locus4->object_id()] = locus4;

    // LocusList specifications

    Parameters parameters_ll_all_A;
    parameters_ll_all_A.insert_name_value("chromosome:position", "5 50000");
    parameters_ll_all_A.insert_name_value("chromosome:position", "6 60000");

    LocusListPtr ll_all_A(new LocusList("ll_all_A"));
    ll_all_A->configure(parameters_ll_all_A, registry);

    Parameters parameters_ll_all_B;
    parameters_ll_all_B.insert_name_value("loci", "locus4");

    LocusListPtr ll_all_B(new LocusList("ll_all_B"));
    ll_all_B->configure(parameters_ll_all_B, registry);

    Parameters parameters_ll_1;
    parameters_ll_1.insert_name_value("loci", "locus1 locus2");

    LocusListPtr ll_1(new LocusList("ll_1"));
    ll_1->configure(parameters_ll_1, registry);

    Parameters parameters_ll_2_A;
    parameters_ll_2_A.insert_name_value("chromosome:position", "9 90000");
    parameters_ll_2_A.insert_name_value("chromosome:position", "10 100000");

    LocusListPtr ll_2_A(new LocusList("ll_2_A"));
    ll_2_A->configure(parameters_ll_2_A, registry);

    Parameters parameters_ll_2_B;
    parameters_ll_2_B.insert_name_value("loci", "locus3");

    LocusListPtr ll_2_B(new LocusList("ll_2_B"));
    ll_2_B->configure(parameters_ll_2_B, registry);

    Parameters parameters_ll_3;
    parameters_ll_3.insert_name_value("chromosome:position", "7 70000");
    parameters_ll_3.insert_name_value("chromosome:position", "8 80000");

    LocusListPtr ll_3(new LocusList("ll_3"));
    ll_3->configure(parameters_ll_3, registry);

    registry[ll_all_A->object_id()] = ll_all_A;
    registry[ll_all_B->object_id()] = ll_all_B;
    registry[ll_1->object_id()] = ll_1;
    registry[ll_2_A->object_id()] = ll_2_A;
    registry[ll_2_B->object_id()] = ll_2_B;
    registry[ll_3->object_id()] = ll_3;

    // VI_Random specification

    Parameters parameters;

    parameters.insert_name_value("locus_list:population:frequencies", "ll_all_A * 0.001 0.0005 ");
    parameters.insert_name_value("locus_list:population:frequencies", "ll_all_B * 0.004 ");
    parameters.insert_name_value("locus_list:population:frequencies", "ll_1 1 0.001 0.002 ");
    parameters.insert_name_value("locus_list:population:frequencies", "ll_2_A 2 0.001 0.0005 ");
    parameters.insert_name_value("locus_list:population:frequencies", "ll_2_B 2 0.003 ");
    parameters.insert_name_value("locus_list:population:frequencies", "ll_3 3 0.001 0.0005 ");

    VariantIndicator_Random vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    shared_ptr<PCG_TestHW> pcg(new PCG_TestHW); // 3 populations: sizes 1000, 2000, 3000
    Population::Configs popconfigs = pcg->population_configs(0, PopulationDataPtrs());

    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    vi.initialize(simconfig);

    if (os_)
    {
        *os_ << "configuration:\n\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    VariantIndicator_Random::EntryMap entries = vi.entries();

    if (os_)
    {
        *os_ << "EntryMap:\n";
        for (VariantIndicator_Random::EntryMap::const_iterator it=entries.begin(); it!=entries.end(); ++it)
        {
            *os_ << it->first.object_id() << ": ";
            copy(it->second.ids.begin(), it->second.ids.end(), ostream_iterator<unsigned int>(*os_, " "));
            *os_ << endl;
        }
    }

    unit_assert(entries.size() == 10);
    unit_assert(entries.at(*locus1).ids.size() == 2);
    unit_assert(entries.at(*locus2).ids.size() == 4);
    unit_assert(entries.at(*locus3).ids.size() == 12);
    unit_assert(entries.at(*locus4).ids.size() == 48);
    unit_assert(entries.at((*ll_all_A)[0]).ids.size() == 12);
    unit_assert(entries.at((*ll_all_A)[1]).ids.size() == 6);
    unit_assert(entries.at((*ll_2_A)[0]).ids.size() == 4);
    unit_assert(entries.at((*ll_2_A)[1]).ids.size() == 2);
    unit_assert(entries.at((*ll_3)[0]).ids.size() == 6);
    unit_assert(entries.at((*ll_3)[1]).ids.size() == 3);

    //vi.write_file("forqs.vi.testing.txt");
}


void test_VariantIndicator_Random_2()
{
    if (os_) *os_ << "test_VariantIndicator_Random_2()\n\n";

    Configurable::Registry registry;

    Parameters parameters_ll_A;
    parameters_ll_A.insert_name_value("chromosome:position", "1 10000");
    parameters_ll_A.insert_name_value("chromosome:position", "2 20000");
    parameters_ll_A.insert_name_value("chromosome:position", "3 30000");

    LocusListPtr ll_A(new LocusList("ll_A"));
    ll_A->configure(parameters_ll_A, registry);

    Parameters parameters_ll_B;
    parameters_ll_B.insert_name_value("chromosome:position", "4 40000");
    parameters_ll_B.insert_name_value("chromosome:position", "5 50000");
    parameters_ll_B.insert_name_value("chromosome:position", "6 60000");

    LocusListPtr ll_B(new LocusList("ll_B"));
    ll_B->configure(parameters_ll_B, registry);

    registry[ll_A->object_id()] = ll_A;
    registry[ll_B->object_id()] = ll_B;

    Random::DistributionPtr distribution = Random::create_constant_distribution("dist", .0005);
    registry[distribution->object_id()] = distribution;

    // VI_Random specification

    Parameters parameters;

    parameters.insert_name_value("locus_list:population:frequencies", "ll_A * 0.001 0.001 0.001 ");
    parameters.insert_name_value("locus_list:population:frequency_distribution", "ll_B * dist");
    
    VariantIndicator_Random vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    shared_ptr<PCG_TestHW> pcg(new PCG_TestHW); // 3 populations: sizes 1000, 2000, 3000
    Population::Configs popconfigs = pcg->population_configs(0, PopulationDataPtrs());

    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    vi.initialize(simconfig);

    if (os_)
    {
        *os_ << "configuration:\n\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    VariantIndicator_Random::EntryMap entries = vi.entries();

    if (os_)
    {
        *os_ << "EntryMap:\n";
        for (VariantIndicator_Random::EntryMap::const_iterator it=entries.begin(); it!=entries.end(); ++it)
        {
            *os_ << it->first.object_id() << ": ";
            copy(it->second.ids.begin(), it->second.ids.end(), ostream_iterator<unsigned int>(*os_, " "));
            *os_ << endl;
        }
    }

    unit_assert(entries.size() == 6);
    unit_assert(entries.at((*ll_A)[0]).ids.size() == 12);
    unit_assert(entries.at((*ll_A)[1]).ids.size() == 12);
    unit_assert(entries.at((*ll_A)[2]).ids.size() == 12);
    unit_assert(entries.at((*ll_B)[0]).ids.size() == 6);
    unit_assert(entries.at((*ll_B)[1]).ids.size() == 6);
    unit_assert(entries.at((*ll_B)[2]).ids.size() == 6);
 
    //vi.write_file("forqs.vi.testing.txt");
}


void test_VariantIndicator_File()
{
    if (os_) *os_ << "test_VariantIndicator_File()\n";

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

    Parameters parameters;
    parameters.insert_name_value("msfile", "../examples/ms_test_data_3.txt");
    parameters.insert_name_value("loci", "locus1 locus_list ");

    VariantIndicator_File vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    for (unsigned int id=0; id<12; ++id)
    {
        if (os_) *os_ << id << " " << vi(id, *locus1) << " " << vi(id, *locus2) << " " << vi(id, *locus3) << endl;

        unit_assert(vi(id, *locus1) == int(id%3==0));
        unit_assert(vi(id, *locus2) == int(id%3==1));
        unit_assert(vi(id, *locus3) == int(id%3==2));
    }
}


void test_VariantIndicator_Mutable()
{
    if (os_) *os_ << "test_VariantIndicator_Mutable()\n";

    // instantiate loci and internal VariantIndicator

    LocusPtr locus1(new Locus("locus1", 0, 100000));
    LocusPtr locus2(new Locus("locus2", 0, 200000));
    LocusPtr locus3(new Locus("locus3", 0, 300000));

    Configurable::Registry registry;
    registry["locus1"] = locus1;
    registry["locus2"] = locus2;
    registry["locus3"] = locus3;

    Parameters parameters_vi_id_range;
    parameters_vi_id_range.insert_name_value("locus:start:count:step:value", "locus1 0 10 1 1");

    VariantIndicatorPtr vi_id_range(new VariantIndicator_IDRange("vi_id_range"));
    vi_id_range->configure(parameters_vi_id_range, registry);

    registry["vi_id_range"] = vi_id_range;

    // instantiate VI_Mutable

    const unsigned int unused_id_start = 1000;

    Parameters parameters;
    parameters.insert_name_value("unused_id_start", unused_id_start);
    parameters.insert_name_value("variant_indicator", "vi_id_range");
    
    VariantIndicator_Mutable vi("id_dummy");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    // test mutate()

    for (unsigned int id=0; id<10; ++id)
        unit_assert(vi(id, *locus1) == 1 && vi(id, *locus2) == 0 && vi(id, *locus3) == 0);

    for (unsigned int id=10; id<20; ++id)
        unit_assert(vi(id, *locus1) == 0 && vi(id, *locus2) == 0 && vi(id, *locus3) == 0);

    unsigned int id_mutant_1 = vi.mutate(0, *locus1, 2);
    if (os_) *os_ << "id_mutant_1: " << id_mutant_1 << endl;
    unit_assert(id_mutant_1 == unused_id_start);
    unit_assert(vi(id_mutant_1, *locus1) == 2);
    unit_assert(vi(id_mutant_1, *locus2) == 0);
    unit_assert(vi(id_mutant_1, *locus3) == 0);

    Loci loci = vi.loci(0, true);
    unit_assert(loci.size() == 1);
    unit_assert(loci.count(*locus1));

    unsigned int id_mutant_2 = vi.mutate(9, *locus2, 3);
    if (os_) *os_ << "id_mutant_2: " << id_mutant_2 << endl;
    unit_assert(vi(id_mutant_2, *locus1) == 1);
    unit_assert(vi(id_mutant_2, *locus2) == 3);
    unit_assert(vi(id_mutant_2, *locus3) == 0);

    loci = vi.loci(0, true);
    unit_assert(loci.size() == 2);
    unit_assert(loci.count(*locus1));
    unit_assert(loci.count(*locus2));

    unsigned int id_mutant_3 = vi.mutate(2, *locus3, 4);
    if (os_) *os_ << "id_mutant_3: " << id_mutant_3 << endl;
    unit_assert(vi(id_mutant_3, *locus1) == 1);
    unit_assert(vi(id_mutant_3, *locus2) == 0);
    unit_assert(vi(id_mutant_3, *locus3) == 4);

    loci = vi.loci(0, true);
    unit_assert(loci.size() == 3);
    unit_assert(loci.count(*locus1));
    unit_assert(loci.count(*locus2));
    unit_assert(loci.count(*locus3));

    if (os_) 
    {
        *os_ << "loci: ";
        copy(loci.begin(), loci.end(), ostream_iterator<Locus>(*os_, " "));
        *os_ << endl;
        vi.report(*os_);
    }

    // make sure nothing has changed with ids 0-19

    for (unsigned int id=0; id<10; ++id)
        unit_assert(vi(id, *locus1) == 1 && vi(id, *locus2) == 0 && vi(id, *locus3) == 0);

    for (unsigned int id=10; id<20; ++id)
        unit_assert(vi(id, *locus1) == 0 && vi(id, *locus2) == 0 && vi(id, *locus3) == 0);
}


class VariantIndicator_Test : public VariantIndicator
{
    public:

    VariantIndicator_Test(const string& id, unsigned int chunk_id, 
                          const Locus& locus, unsigned int value) 
    :   Configurable(id), chunk_id_(chunk_id), locus_(locus), value_(value)
    {} 

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const 
    {
        return chunk_id == chunk_id_ && locus_ == locus ? value_ : 0;
    }

    virtual std::string class_name() const {return "VariantIndicator_Test";}
    virtual Parameters parameters() const {return Parameters();}

    private:

    unsigned int chunk_id_;
    Locus locus_;
    unsigned int value_;
};


void test_VariantIndicator_Composite()
{
    Locus locus1("locus1", 1, 1000);
    Locus locus2("locus2", 2, 2000);
    Locus locus3("locus3", 3, 3000);

    VariantIndicatorPtr vi1(new VariantIndicator_Test("vi1", 10, locus1, 1));
    VariantIndicatorPtr vi2(new VariantIndicator_Test("vi2", 20, locus2, 2));
    VariantIndicatorPtr vi3(new VariantIndicator_Test("vi3", 30, locus3, 3));

    Configurable::Registry registry;
    registry["vi1"] = vi1;
    registry["vi2"] = vi2;
    registry["vi3"] = vi3;

    Parameters parameters;
    parameters.insert_name_value("variant_indicators", "vi1 vi2 vi3 ");

    VariantIndicator_Composite vi("vi");
    vi.configure(parameters, registry);

    Parameters parameters_out = vi.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "parameters_out:\n" << parameters_out << endl;

        *os_ << "configuration:\n";
        set<string> ids_written;
        vi.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    unit_assert(vi(10, locus1) == 1);
    unit_assert(vi(11, locus1) == 0);
    unit_assert(vi(10, locus2) == 0);
    unit_assert(vi(10, locus3) == 0);

    unit_assert(vi(20, locus2) == 2);
    unit_assert(vi(21, locus2) == 0);
    unit_assert(vi(20, locus1) == 0);
    unit_assert(vi(20, locus3) == 0);

    unit_assert(vi(30, locus3) == 3);
    unit_assert(vi(31, locus3) == 0);
    unit_assert(vi(30, locus1) == 0);
    unit_assert(vi(30, locus2) == 0);
}


void test()
{
    test_VariantIndicator_SingleLocusHardyWeinberg();
    test_Configurable_VariantIndicator_SingleLocusHardyWeinberg();
    test_VariantIndicator_TwoLocusLD();
    test_VariantIndicator_IDRange();
    test_VariantIndicator_IDSet();
    test_VariantIndicator_Random();
    test_VariantIndicator_Random_2();
    test_VariantIndicator_File();
    test_VariantIndicator_Mutable();
    test_VariantIndicator_Composite();
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


