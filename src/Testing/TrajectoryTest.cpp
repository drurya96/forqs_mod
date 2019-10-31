//
// TrajectoryTest.cpp
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


#include "Trajectory.hpp"
#include "PopulationConfigGenerator.hpp"
#include "Simulator.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


class Trajectory_TestingComposite : public Trajectory
{
    public:

    Trajectory_TestingComposite(const string& id) : Trajectory(id) {}

    virtual double value(size_t generation_index,
                         size_t population_index) const
    {
        population_indices.push_back(population_index);        
        return 0;
    }

    virtual std::string class_name() const {return "Trajectory_TestingComposite";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}

    mutable vector<size_t> population_indices;
};


typedef boost::shared_ptr<Trajectory_TestingComposite> Trajectory_TestingCompositePtr;


class PopulationConfigGenerator_Dummy : public PopulationConfigGenerator
{
    public:

    PopulationConfigGenerator_Dummy(size_t population_count)
    :   PopulationConfigGenerator("dummy")
    {
        population_count_ = population_count; // protected
    }

    virtual Population::Configs population_configs(size_t generation_index,
                                                   const PopulationDataPtrs&) const
    {
        return Population::Configs(population_count_);
    }
};

 
void test_Trajectory_PopulationComposite()
{
    if (os_) *os_ << "test_Trajectory_PopulationComposite()\n";

    string id_testing_0 = "id_testing_0";
    string id_testing_1 = "id_testing_1";
    string id_testing_2 = "id_testing_2";

    TrajectoryPtr trajectory_testing_0(new Trajectory_TestingComposite(id_testing_0));
    TrajectoryPtr trajectory_testing_1(new Trajectory_TestingComposite(id_testing_1));
    TrajectoryPtr trajectory_testing_2(new Trajectory_TestingComposite(id_testing_2));

    ostringstream trajectories_string;
    trajectories_string << id_testing_0 << " " << id_testing_1 << " " << id_testing_2 << " ";

    Parameters parameters;
    parameters.insert_name_value("trajectories", trajectories_string.str());

    Configurable::Registry registry;
    registry[id_testing_0] = trajectory_testing_0;
    registry[id_testing_1] = trajectory_testing_1;
    registry[id_testing_2] = trajectory_testing_2;

    Trajectory_PopulationComposite trajectory("trajectory_population_composite");
    trajectory.configure(parameters, registry);

    Parameters parameters_out = trajectory.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
        
        *os_ << "configuration:\n";
        set<string> ids_written;
        trajectory.write_configuration(*os_, ids_written);
    }

    unit_assert(parameters == parameters_out);

    PopulationConfigGeneratorPtr pcg(new PopulationConfigGenerator_Dummy(3)); // 3 populations
    SimulatorConfig simconfig;
    simconfig.population_config_generator = pcg;

    trajectory.initialize(simconfig);

    trajectory.value(0, 0);
    trajectory.value(0, 1);
    trajectory.value(0, 1);
    trajectory.value(0, 2);
    trajectory.value(0, 2);
    trajectory.value(0, 2);

    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_0)
            ->population_indices.size() == 1);
    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_1)
            ->population_indices.size() == 2);
    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_2)
            ->population_indices.size() == 3);
}


void test_Trajectory_GenerationComposite()
{
    if (os_) *os_ << "test_Trajectory_GenerationComposite()\n";

    string id_testing_0 = "id_testing_0";
    string id_testing_1 = "id_testing_1";
    string id_testing_2 = "id_testing_2";

    TrajectoryPtr trajectory_testing_0(new Trajectory_TestingComposite(id_testing_0));
    TrajectoryPtr trajectory_testing_1(new Trajectory_TestingComposite(id_testing_1));
    TrajectoryPtr trajectory_testing_2(new Trajectory_TestingComposite(id_testing_2));

    Parameters parameters;
    parameters.insert_name_value("generation:trajectory", "0 " + id_testing_0);
    parameters.insert_name_value("generation:trajectory", "5 " + id_testing_1);
    parameters.insert_name_value("generation:trajectory", "23 " + id_testing_2);

    Configurable::Registry registry;
    registry[id_testing_0] = trajectory_testing_0;
    registry[id_testing_1] = trajectory_testing_1;
    registry[id_testing_2] = trajectory_testing_2;

    Trajectory_GenerationComposite trajectory("trajectory_population_composite");
    trajectory.configure(parameters, registry);

    Parameters parameters_out = trajectory.parameters();

    set<string> ids_written;
    if (os_) trajectory.write_configuration(*os_, ids_written);

    unit_assert(parameters == parameters_out);

    for (size_t generation_index=0; generation_index<30; ++generation_index)
        trajectory.value(generation_index, 0);

    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_0)
            ->population_indices.size() == 5);
    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_1)
            ->population_indices.size() == 18);
    unit_assert(dynamic_pointer_cast<Trajectory_TestingComposite>(trajectory_testing_2)
            ->population_indices.size() == 7);
}


void test_Trajectory_Constant()
{
    if (os_) *os_ << "test_Trajectory_Constant()\n";

    const double value_0 = 1.23;
    const double value_1 = 2.34;
    const double value_2 = 3.45;
    const double values[] = {value_0, value_1, value_2};

    Parameters parameters_0;
    parameters_0.insert_name_value("value", value_0);

    Configurable::Registry registry;
    TrajectoryPtr t0(new Trajectory_Constant("id0"));
    t0->configure(parameters_0, registry);

    unit_assert(t0->value(0,0) == value_0);

    Parameters parameters_0_out = t0->parameters();
    unit_assert(parameters_0 == parameters_0_out);

    TrajectoryPtr t1(new Trajectory_Constant("id1", value_1));
    TrajectoryPtr t2(new Trajectory_Constant("id2", value_2));
     
    // additional testing with PopulationComposite

    TrajectoryPtrs ts;
    ts.push_back(t0);
    ts.push_back(t1);
    ts.push_back(t2);

    TrajectoryPtr t_012(new Trajectory_PopulationComposite("id012", ts));

    set<string> written;
    if (os_) t_012->write_configuration(*os_, written);

    for (size_t generation_index=0; generation_index<20; ++generation_index)
    {
        for (size_t population_index=0; population_index<3; ++population_index)
        {
            if (os_) *os_ << t_012->value(generation_index, population_index) << " ";
            unit_assert(t_012->value(generation_index, population_index) == values[population_index]);
        }
        if (os_) *os_ << endl;
    }
}


void test_Trajectory_Linear()
{
    const double m = 2.0;
    const double b = 1.0;
    
    TrajectoryPtr t(new Trajectory_Linear("id_linear", m, b ));
    
    unit_assert(t->value(0, 234) == b);
    unit_assert(t->value(1, 234) == b + m);
    unit_assert(t->value(2, 234) == b + 2*m);

    // test alternate parameterization

    Parameters parameters;
    parameters.insert_name_value("begin:value", "5 0");
    parameters.insert_name_value("end:value", "10 50");

    Configurable::Registry registry;
    t->configure(parameters, registry);

    unit_assert(t->value(5, 456) == 0);
    unit_assert(t->value(6, 456) == 10);
    unit_assert(t->value(7, 456) == 20);
    unit_assert(t->value(8, 456) == 30);
    unit_assert(t->value(9, 456) == 40);
    unit_assert(t->value(10, 456) == 50);

    set<string> written;
    if (os_) t->write_configuration(*os_, written);

    Parameters parameters_out = t->parameters();
    unit_assert(parameters == parameters_out);
}


void test_Trajectory_Exponential()
{
    if (os_) *os_ << "test_Trajectory_Exponential()\n";

    const double t0 = 5;
    const double A = 32;
    const double r = log(2);

    Parameters parameters;
    parameters.insert_name_value("generation_begin", t0);
    parameters.insert_name_value("value_begin", A);
    parameters.insert_name_value("rate", r);
    
    Trajectory_Exponential trajectory("id_dummy");
    Configurable::Registry registry;
    trajectory.configure(parameters, registry);

    Parameters parameters_out = trajectory.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
        *os_ << "configuration:\n";
        set<string> written;
        trajectory.write_configuration(*os_, written);

        for (size_t generation_index=5; generation_index<10; ++generation_index)
            *os_ << generation_index << " " << trajectory.value(generation_index, 0) << endl;
    }

    unit_assert(parameters == parameters_out);

    const double epsilon = 1e-6;
    unit_assert_equal(trajectory.value(5,0), 32, epsilon);
    unit_assert_equal(trajectory.value(6,1), 64, epsilon);
    unit_assert_equal(trajectory.value(7,1000), 128, epsilon);
    unit_assert_equal(trajectory.value(8,0), 256, epsilon);
    unit_assert_equal(trajectory.value(9,0), 512, epsilon);
    unit_assert_equal(trajectory.value(10,0), 1024, epsilon);
}


void test()
{
    test_Trajectory_PopulationComposite();
    test_Trajectory_GenerationComposite();
    test_Trajectory_Constant();
    test_Trajectory_Linear();
    test_Trajectory_Exponential();
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


