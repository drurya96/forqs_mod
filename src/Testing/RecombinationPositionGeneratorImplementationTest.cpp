//
// RecombinationPositionGeneratorImplementationTest.cpp
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


#include "RecombinationPositionGeneratorImplementation.hpp"
#include "PopulationConfigGenerator.hpp"
#include "Simulator.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_RecombinationPositionGenerator_Trivial()
{
    if (os_) *os_ << "test_RecombinationPositionGenerator_Trivial()\n";

    RecombinationPositionGenerator_Trivial r("dummy_id");

    size_t count_total = 10000;
    double count_empty = 0;
    for (size_t i=0; i<count_total; i++)
    {
        vector<unsigned int> positions = r.get_positions(0);
        count_empty += positions.empty(); 
    }
    
    double frequency = count_empty / count_total;
    double epsilon = 3. / 2 / 100; // 3 standard deviations

    if (os_) *os_ << "frequency: " << frequency << endl
                  << "epsilon: " << epsilon << endl;

    unit_assert_equal(frequency, .5, epsilon);
}


class PopulationConfigGenerator_ChromosomeLengths : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_ChromosomeLengths()
    :   PopulationConfigGenerator("dummy")
    {
        chromosome_lengths_.push_back(1000);
        chromosome_lengths_.push_back(20000);
        chromosome_lengths_.push_back(300000);
        chromosome_pair_count_ = chromosome_lengths_.size();
    }

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs& population_datas) const
    {
        return Population::Configs();
    }
};


void test_RecombinationPositionGenerator_SingleCrossover()
{
    if (os_) *os_ << "test_RecombinationPositionGenerator_SingleCrossover()\n";
    RecombinationPositionGenerator_SingleCrossover uniform("dummy_id");
    RecombinationPositionGenerator& r = uniform; // default param in base interface only

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_ChromosomeLengths);
    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    r.initialize(simconfig);

    size_t count_total = 1000;

    double count1 = 0;
    double count1_recombined = 0;
    double count2 = 0;
    double count2_recombined = 0;

    for (size_t i=0; i<count_total; i++)
    {
        vector<unsigned int> positions = r.get_positions();
        
        if (positions.empty())
            ++count1;
        else if (positions.size() == 1)
        {
            if (positions[0] == 0)
                ++count2;
            else
                ++count1_recombined;
        }
        else
        {
            unit_assert(positions.size() == 2 && positions[0] == 0); 
            ++count2_recombined;
        }
    }

    vector<double> v;
    v.push_back(count1);
    v.push_back(count1_recombined);
    v.push_back(count2);
    v.push_back(count2_recombined);

    double chisq = 0;
    double e = count_total/4.;
    for (vector<double>::const_iterator it=v.begin(); it!=v.end(); ++it)
        chisq += (*it - e) * (*it - e) / e;

    if (os_) *os_ << "count1: " << count1 << endl
                  << "count1_recombined: " << count1_recombined << endl
                  << "count2: " << count2 << endl
                  << "count2_recombined: " << count2_recombined << endl
                  << "chisq: " << chisq << endl;
    
    double cutoff = 11.34487; // from R: qchisq(.99, df=3)

    unit_assert(chisq < cutoff);
}


void test_Configurable_RecombinationPositionGenerator_RecombinationMap()
{
    Parameters parameters_in;
    parameters_in.insert_name_value("filename", "../examples/genetic_map_chr21_b36.txt");
    parameters_in.insert_name_value("filename", "../examples/genetic_map_chr21_b36.txt");
    parameters_in.insert_name_value("filename", "../examples/genetic_map_chr21_b36.txt");

    RecombinationPositionGenerator_RecombinationMap rpg("dummy_id");

    Configurable::Registry registry;
    rpg.configure(parameters_in, registry);

    Parameters parameters_out = rpg.parameters();
    unit_assert(parameters_in == parameters_out);

    vector<string> filenames = parameters_out.values<string>("filename");
    unit_assert(filenames.size() == 3);
    for (size_t i=0; i<3; ++i)
        unit_assert(filenames[i] == "../examples/genetic_map_chr21_b36.txt");

    if (os_)
    {
        *os_ << "filenames: ";
        copy(filenames.begin(), filenames.end(), ostream_iterator<string>(*os_, " "));
        *os_ << endl;
    }
}


string positions_to_string(const vector<unsigned int>& positions)
{
    ostringstream oss;
    oss << "(";
    copy(positions.begin(), positions.end(), ostream_iterator<unsigned int>(oss, " "));
    oss << ") ";
    return oss.str();
}


void test_Configurable_RecombinationPositionGenerator_Uniform()
{
    if (os_) *os_ << "test_Configurable_RecombinationPositionGenerator_Uniform()\n";

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_ChromosomeLengths);
    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;
    
    Parameters parameters_in;
    parameters_in.insert_name_value("rates", "1 2 3 ");

    RecombinationPositionGenerator_Uniform rpg("dummy_id");
    Configurable::Registry registry;
    rpg.configure(parameters_in, registry);
    rpg.initialize(simconfig);

    Parameters parameters_out = rpg.parameters();

    if (os_)
    {
        *os_ << "parameters_in:\n" << parameters_in << endl
             << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters_in == parameters_out);

    // demo

    if (os_)
    {
        size_t maternal_count = 0, paternal_count = 0;
        vector<double> recombination_event_counts(3);
        const size_t replicate_count = 10;

        for (size_t i=0; i<replicate_count; ++i)
        {
            for (int chromosome_pair_index=0; chromosome_pair_index<3; ++chromosome_pair_index)
            {
                vector<unsigned int> positions = rpg.get_positions(chromosome_pair_index);
                *os_ << positions_to_string(positions);

                size_t position_count = positions.size();

                if (!positions.empty() && positions[0]==0)
                {
                    --position_count;
                    ++maternal_count;
                }
                else
                {
                    ++paternal_count;
                }

                recombination_event_counts[chromosome_pair_index] += position_count;
            }
            *os_ << endl;
        }

        *os_ << "maternal_count: " << maternal_count << endl;
        *os_ << "paternal_count: " << paternal_count << endl;
        *os_ << "average # of recombination events: ";
        for (size_t chromosome_pair_index=0; chromosome_pair_index<3; ++chromosome_pair_index)
            *os_ << recombination_event_counts[chromosome_pair_index]/replicate_count << " ";
        *os_ << endl;
    }
}


class RPGCounter : public RecombinationPositionGenerator
{ 
    public:

    RPGCounter(const std::string& id)
    :   RecombinationPositionGenerator(id)
    {}

    virtual std::vector<unsigned int> get_positions(size_t chromosome_pair_index) const
    {
        ++counter[chromosome_pair_index];
        return vector<unsigned int>();
    }

    //virtual std::string class_name() const {return "RecombinationPositionGenerator_Counter";}

    typedef map<size_t,unsigned int> Counter;
    mutable Counter counter;
};


void test_RecombinationPositionGenerator_Composite()
{
    boost::shared_ptr<RPGCounter> rpg_default(new RPGCounter("rpg_default"));
    boost::shared_ptr<RPGCounter> rpg_3(new RPGCounter("rpg_3"));
    boost::shared_ptr<RPGCounter> rpg_5(new RPGCounter("rpg_5"));

    Configurable::Registry registry;
    registry["rpg_default"] = rpg_default;
    registry["rpg_3"] = rpg_3;
    registry["rpg_5"] = rpg_5;

    Parameters parameters_in;
    parameters_in.insert_name_value("default_recombination_position_generator", "rpg_default");
    parameters_in.insert_name_value("chromosome:recombination_position_generator", "3 rpg_3");
    parameters_in.insert_name_value("chromosome:recombination_position_generator", "5 rpg_5");

    RecombinationPositionGenerator_Composite rpg("rpg");
    rpg.configure(parameters_in, registry);

    Parameters parameters_out = rpg.parameters();

    if (os_)
    {
        *os_ << "parameters in:\n" << parameters_in << endl;
        *os_ << "parameters out:\n" << parameters_out << endl;
    }

    unit_assert(parameters_in == parameters_out);

    for (size_t i=0; i<10; ++i) rpg.get_positions(i);
    for (size_t i=0; i<10; ++i) rpg.get_positions(i);

    unit_assert(rpg_3->counter.size() == 1 && rpg_3->counter[2] == 2);
    unit_assert(rpg_5->counter.size() == 1 && rpg_5->counter[4] == 2);
    unit_assert(rpg_default->counter.size() == 8);
    for (size_t i=0; i<10; ++i)
        if (i!=2 && i!=4)
            unit_assert(rpg_default->counter[i] == 2);
}


void test()
{
    test_RecombinationPositionGenerator_Trivial();
    test_RecombinationPositionGenerator_SingleCrossover();
    test_Configurable_RecombinationPositionGenerator_RecombinationMap();
    test_Configurable_RecombinationPositionGenerator_Uniform();
    test_RecombinationPositionGenerator_Composite();
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


