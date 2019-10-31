//
// LocusTest.cpp
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


#include "Locus.hpp"
#include "PopulationConfigGenerator.hpp"
#include "Simulator.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_Configurable_Locus()
{
    if (os_) *os_ << "test_Configurable_Locus()\n";

    Parameters parameters_in;
    parameters_in.insert_name_value("chromosome", 1);
    parameters_in.insert_name_value("position", 123456);

    if (os_) *os_ << "parameters_in:\n" << parameters_in << endl;

    Locus locus("my_locus");

    Configurable::Registry registry;
    locus.configure(parameters_in, registry);

    Parameters parameters_out = locus.parameters();

    if (os_) *os_ << "parameters_out:\n" << parameters_out << endl;

    unit_assert(parameters_in == parameters_out);
}


void test_Configurable_LocusList()
{
    if (os_) *os_ << "test_Configurable_LocusList()\n";

    Configurable::Registry registry;
    LocusPtr locus(new Locus("my_locus", 2, 234567));
    LocusPtr locus2(new Locus("my_locus_2", 4, 456789));
    LocusPtr locus3(new Locus("my_locus_3", 6, 678901));
    LocusPtr locus4(new Locus("my_locus_4", 8, 890123));
    registry[locus->object_id()] = locus;
    registry[locus2->object_id()] = locus2;
    registry[locus3->object_id()] = locus3;
    registry[locus4->object_id()] = locus4;

    LocusList locus_list("locus_list");

    // unnamed "chromosome:position"

    Parameters parameters_chr;
    parameters_chr.insert_name_value("chromosome:position", "1 123456");
    parameters_chr.insert_name_value("chromosome:position", "3 345678");

    locus_list.configure(parameters_chr, registry);

    Parameters parameters_out = locus_list.parameters();

    if (os_) 
    {
        *os_ << "parameters_chr:\n" << parameters_chr << endl;

        *os_ << "configuration:\n";
        set<string> written;
        locus_list.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters_chr == parameters_out);
    unit_assert(locus_list.size() == 2);

    // named "loci"

    Parameters parameters_loci;
    parameters_loci.insert_name_value("loci", "my_locus my_locus_2 ");
    parameters_loci.insert_name_value("loci", "my_locus_3 my_locus_4 ");

    locus_list.configure(parameters_loci, registry);

    parameters_out = locus_list.parameters();

    if (os_) 
    {
        *os_ << "parameters_loci:\n" << parameters_loci << endl;

        *os_ << "configuration:\n";
        set<string> written;
        locus_list.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(locus_list.size() == 4);

    // specification of both "chromosome:position" and "loci":
    //  all unnamed loci are placed first in the array

    LocusListPtr sublist(new LocusList("sublist"));
    sublist->push_back(*locus3);
    sublist->push_back(*locus4);
    registry[sublist->object_id()] = sublist;

    Parameters parameters_both;
    parameters_both.insert_name_value("chromosome:position", "1 123456");
    parameters_both.insert_name_value("loci", "my_locus my_locus_2");
    parameters_both.insert_name_value("chromosome:position", "3 345678");
    parameters_both.insert_name_value("loci", "sublist");

    locus_list.configure(parameters_both, registry);

    parameters_out = locus_list.parameters();

    if (os_) 
    {
        *os_ << "parameters_both:\n" << parameters_both << endl;

        *os_ << "configuration:\n";
        set<string> written;
        locus_list.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(locus_list.size() == 6);
    unit_assert(locus_list[0].object_id() == "locus_list[0]");
    unit_assert(locus_list[1].object_id() == "locus_list[1]");
    unit_assert(locus_list[2].object_id() == "my_locus");
    unit_assert(locus_list[3].object_id() == "my_locus_2");
    unit_assert(locus_list[4].object_id() == "my_locus_3");
    unit_assert(locus_list[5].object_id() == "my_locus_4");
}


class PopulationConfigGenerator_ChromosomeLengths : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_ChromosomeLengths()
    :   PopulationConfigGenerator("dummy")
    {
        chromosome_lengths_.push_back(10000);
        chromosome_lengths_.push_back(20000);
        chromosome_lengths_.push_back(30000);
        chromosome_pair_count_ = chromosome_lengths_.size();
    }

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const
    {
        return Population::Configs();
    }
};


void test_Configurable_LocusList_Random()
{
    if (os_) *os_ << "test_Configurable_LocusList_Random()\n";

    Parameters parameters_in;
    parameters_in.insert_name_value("locus_count", 600);

    Configurable::Registry registry;

    LocusList_Random locus_list("locus_list");
    locus_list.configure(parameters_in, registry);

    Parameters parameters_out = locus_list.parameters();

    if (os_) 
    {
        *os_ << "parameters_in:\n" << parameters_in << endl;

        *os_ << "configuration:\n";
        set<string> written;
        locus_list.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters_in == parameters_out);

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_ChromosomeLengths);
    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    unit_assert(locus_list.empty());
    locus_list.initialize(simconfig);

    // if (os_) copy(locus_list.begin(), locus_list.end(), ostream_iterator<Locus>(*os_, "\n"));

    unit_assert(locus_list.size() == 600);

    vector<unsigned int> counts(3);
    for (LocusList::const_iterator it=locus_list.begin(); it!=locus_list.end(); ++it)
        ++counts[it->chromosome_pair_index];

    if (os_)
    {
        *os_ << "counts: ";
        copy(counts.begin(), counts.end(), ostream_iterator<unsigned int>(*os_, " "));
        *os_ << endl;
    }
}


void test()
{
    test_Configurable_Locus();
    test_Configurable_LocusList();
    test_Configurable_LocusList_Random();
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


