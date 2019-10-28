//
// FitnessFunctionImplementationTest.cpp
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


#include "FitnessFunctionImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_FitnessFunction_Optimum_Polynomial()
{
    if (os_) *os_ << "test_FitnessFunction_Optimum_Polynomial()\n";

    // check configuration

    string qtid = "qtid";

    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid);
    parameters.insert_name_value("optimum", 100);
    parameters.insert_name_value("radius:power", "10 2");

    Configurable::Registry registry;

    string id_ff_optimum("ff_optimum");
    FitnessFunction_Optimum ff(id_ff_optimum);

    ff.configure(parameters, registry);
    
    Parameters parameters_out = ff.parameters();

    unit_assert(parameters == parameters_out);

    // check fitness values

    DataVectorPtr trait_values(new DataVector);
    trait_values->push_back(100); // 1
    trait_values->push_back(80);  // 0
    trait_values->push_back(90);  // 0
    trait_values->push_back(120); // 0
    trait_values->push_back(110); // 0
    trait_values->push_back(105); // (.5)^2
    trait_values->push_back(95);  // (.5)^2
    trait_values->push_back(101); // (.9)^2
    trait_values->push_back(99);  // (.9)^2
    trait_values->push_back(109); // (.1)^2
    trait_values->push_back(91);  // (.1)^2

    PopulationData population_data;
    population_data.generation_index = 0;
    population_data.population_index = 0;
    (*population_data.trait_values)[qtid] = trait_values;

    ff.calculate_trait_values(population_data); 
    DataVectorPtr fitnesses = population_data.trait_values->at(id_ff_optimum);

    if (os_) *os_ << "fitnesses: " << *fitnesses << endl;

    const double epsilon = 1e-6;
    unit_assert(fitnesses->size() == trait_values->size());
    unit_assert(fitnesses->at(0) == 1);
    unit_assert(fitnesses->at(1) == 0);
    unit_assert(fitnesses->at(2) == 0);
    unit_assert(fitnesses->at(3) == 0);
    unit_assert(fitnesses->at(4) == 0);
    unit_assert_equal(fitnesses->at(5), .25, epsilon);
    unit_assert_equal(fitnesses->at(6), .25, epsilon);
    unit_assert_equal(fitnesses->at(7), .81, epsilon);
    unit_assert_equal(fitnesses->at(8), .81, epsilon);
    unit_assert_equal(fitnesses->at(9), .01, epsilon);
    unit_assert_equal(fitnesses->at(10), .01, epsilon);
}


void test_FitnessFunction_Optimum_Gaussian()
{
    if (os_) *os_ << "test_FitnessFunction_Optimum_Gaussian()\n";

    // check configuration

    string qtid = "qtid";

    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid);
    parameters.insert_name_value("optimum", 100);
    parameters.insert_name_value("gaussian_width", "5");

    Configurable::Registry registry;

    string id_ff_optimum("ff_optimum");
    FitnessFunction_Optimum ff(id_ff_optimum);

    ff.configure(parameters, registry);
    
    Parameters parameters_out = ff.parameters();

    unit_assert(parameters == parameters_out);

    // check fitness values

    DataVectorPtr trait_values(new DataVector);
    trait_values->push_back(100); // 1
    trait_values->push_back(105); // exp(-1/2)
    trait_values->push_back(95);  // exp(-1/2)
    trait_values->push_back(110); // exp(-4/2)
    trait_values->push_back(90);  // exp(-4/2)

    PopulationData population_data;
    population_data.generation_index = 0;
    population_data.population_index = 0;
    (*population_data.trait_values)[qtid] = trait_values;

    ff.calculate_trait_values(population_data); 
    DataVectorPtr fitnesses = population_data.trait_values->at(id_ff_optimum);

    if (os_) *os_ << "fitnesses: " << *fitnesses << endl;

    const double epsilon = 1e-6;
    unit_assert(fitnesses->size() == trait_values->size());
    unit_assert(fitnesses->at(0) == 1);
    unit_assert_equal(fitnesses->at(1), exp(-.5), epsilon);
    unit_assert_equal(fitnesses->at(2), exp(-.5), epsilon);
    unit_assert_equal(fitnesses->at(3), exp(-2), epsilon);
    unit_assert_equal(fitnesses->at(4), exp(-2), epsilon);
}


void test_FitnessFunction_TruncationSelection()
{
    if (os_) *os_ << "test_FitnessFunction_TruncationSelection()\n";

    // check configuration

    string qtid = "qtid";

    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid);
    parameters.insert_name_value("proportion_selected", .2);

    Configurable::Registry registry;
    string id_ff("ff");
    FitnessFunction_TruncationSelection ff_truncation(id_ff);
    QuantitativeTrait& ff(ff_truncation);
    ff.configure(parameters, registry);
    
    Parameters parameters_out = ff.parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;

        *os_ << "configuration:\n";
        set<string> written;
        ff.write_configuration(*os_, written);

        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    // check fitness values

    DataVectorPtr trait_values(new DataVector);
    trait_values->push_back(100); // 1
    trait_values->push_back(80);  // 0
    trait_values->push_back(90);  // 0
    trait_values->push_back(120); // 0
    trait_values->push_back(110); // 0
    trait_values->push_back(105); // (.5)^2
    trait_values->push_back(95);  // (.5)^2
    trait_values->push_back(101); // (.9)^2
    trait_values->push_back(99);  // (.9)^2
    trait_values->push_back(109); // (.1)^2
    trait_values->push_back(91);  // (.1)^2


    PopulationDataPtr population_data(new PopulationData);
    population_data->generation_index = 0;
    population_data->population_index = 0;
    (*population_data->trait_values)[qtid] = trait_values;

    PopulationDataPtrs population_datas;
    population_datas.push_back(population_data);

    ff.calculate_trait_values(population_datas); 
    DataVectorPtr fitnesses = population_data->trait_values->at(id_ff);

    if (os_) *os_ << "trait_values: " << *trait_values << endl
                  << "fitnesses: " << *fitnesses << endl;

    unit_assert(fitnesses->size() == trait_values->size());
    for (size_t i=0; i<11; ++i)
        unit_assert(fitnesses->at(i) == (i==3 || i==4 ? 1 : 0));
}


void test_FitnessFunction_TruncationSelection_2()
{
    if (os_) *os_ << "test_FitnessFunction_TruncationSelection_2()\n";

    DataVectorPtr trait_values_1(new DataVector);
    for (int i=0; i<10; ++i) trait_values_1->push_back(i);
    DataVectorPtr trait_values_2(new DataVector);
    for (int i=0; i<10; ++i) trait_values_2->push_back(2*i);

    if (os_)
    {
        *os_ << "trait_values_1: " << *trait_values_1 << endl;
        *os_ << "trait_values_2: " << *trait_values_2 << endl;
    }

    string qtid = "qtid";

    PopulationDataPtr population_data_1(new PopulationData);
    (*population_data_1->trait_values)[qtid] = trait_values_1;

    PopulationDataPtr population_data_2(new PopulationData);
    (*population_data_2->trait_values)[qtid] = trait_values_2;

    PopulationDataPtrs population_datas;
    population_datas.push_back(population_data_1);
    population_datas.push_back(population_data_2);

    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid);
    parameters.insert_name_value("proportion_selected", .5);
    parameters.insert_name_value("lower_tail", 1);

    Configurable::Registry registry;

    FitnessFunction_TruncationSelection ff("ff");
    ff.configure(parameters, registry);

    Parameters parameters_out = ff.parameters();

    if (os_)
    {
        *os_ << "parameters in:\n" << parameters << endl;
        *os_ << "parameters out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    // test multiple thresholds

    ff.calculate_trait_values(population_datas);
    DataVectorPtr ff_values_1 = (*population_datas.front()->trait_values)["ff"];
    DataVectorPtr ff_values_2 = (*population_datas.back()->trait_values)["ff"];

    if (os_)
    {
        *os_ << "multiple thresholds:\n";
        *os_ << "ff_values_1: " << *ff_values_1 << endl;
        *os_ << "ff_values_2: " << *ff_values_2 << endl << endl;
    }

    unit_assert(ff_values_1->size() == 10);
    unit_assert(ff_values_2->size() == 10);

    for (int i=0; i<5; ++i) unit_assert(ff_values_1->at(i) == 1);
    for (int i=5; i<10; ++i) unit_assert(ff_values_1->at(i) == 0);
    for (int i=0; i<5; ++i) unit_assert(ff_values_2->at(i) == 1);
    for (int i=5; i<10; ++i) unit_assert(ff_values_2->at(i) == 0);

    // test single threshold

    population_datas.front()->trait_values->erase("ff");
    population_datas.back()->trait_values->erase("ff");

    parameters.insert_name_value("single_threshold_population", 1);

    ff.configure(parameters, registry);

    ff.calculate_trait_values(population_datas);
    ff_values_1 = (*population_datas.front()->trait_values)["ff"];
    ff_values_2 = (*population_datas.back()->trait_values)["ff"];

    if (os_)
    {
        *os_ << "single threshold:\n";
        *os_ << "ff_values_1: " << *ff_values_1 << endl;
        *os_ << "ff_values_2: " << *ff_values_2 << endl << endl;
    }

    unit_assert(ff_values_1->size() == 10);
    unit_assert(ff_values_2->size() == 10);

    for (int i=0; i<5; ++i) unit_assert(ff_values_1->at(i) == 1);
    for (int i=5; i<10; ++i) unit_assert(ff_values_1->at(i) == 0);
    for (int i=0; i<3; ++i) unit_assert(ff_values_2->at(i) == 1);
    for (int i=3; i<10; ++i) unit_assert(ff_values_2->at(i) == 0);
}


void test_FitnessFunction_TruncationSelection_3()
{
    if (os_) *os_ << "test_FitnessFunction_TruncationSelection_3()\n";

    DataVectorPtr trait_values(new DataVector);
    for (int i=0; i<10; ++i) trait_values->push_back(i+1);
    for (int i=0; i<10; ++i) trait_values->push_back(0);

    if (os_)
        *os_ << "trait_values: " << *trait_values<< endl;

    string qtid = "qtid";

    PopulationDataPtr population_data(new PopulationData);
    (*population_data->trait_values)[qtid] = trait_values;

    PopulationDataPtrs population_datas;
    population_datas.push_back(population_data);

    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid);
    parameters.insert_name_value("proportion_selected", .5);

    Configurable::Registry registry;

    FitnessFunction_TruncationSelection ff("ff");
    ff.configure(parameters, registry);

    Parameters parameters_out = ff.parameters();

    if (os_)
    {
        *os_ << endl
             << "ignore_zero = 0:\n"
             << "parameters in:\n" << parameters << endl
             << "parameters out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    // test ignore_zero = 0

    ff.calculate_trait_values(population_datas);
    DataVectorPtr ff_values = (*population_datas.front()->trait_values)["ff"];

    if (os_)
        *os_ << "ff_values: " << *ff_values << endl;

    unit_assert(ff_values->size() == 20);

    for (int i=0; i<10; ++i) unit_assert(ff_values->at(i) == 1);
    for (int i=10; i<20; ++i) unit_assert(ff_values->at(i) == 0);

    // test ignore_zero = 1

    population_datas.front()->trait_values->erase("ff");
    parameters.insert_name_value("ignore_zero", 1);

    ff.configure(parameters, registry);

    parameters_out = ff.parameters();

    if (os_)
    {
        *os_ << endl 
             << "ignore_zero = 1:\n"
             << "parameters in:\n" << parameters << endl
             << "parameters out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    ff.calculate_trait_values(population_datas);
    ff_values = (*population_datas.front()->trait_values)["ff"];

    if (os_)
        *os_ << "ff_values: " << *ff_values << endl;

    unit_assert(ff_values->size() == 20);

    for (int i=0; i<5; ++i) unit_assert(ff_values->at(i) == 0);
    for (int i=5; i<10; ++i) unit_assert(ff_values->at(i) == 1);
    for (int i=10; i<20; ++i) unit_assert(ff_values->at(i) == 0);

    // test ignore_zero = 1 with lower_tail = 1

    population_datas.front()->trait_values->erase("ff");
    parameters.insert_name_value("ignore_zero", 1);
    parameters.insert_name_value("lower_tail", 1);

    ff.configure(parameters, registry);

    ff.calculate_trait_values(population_datas);
    ff_values = (*population_datas.front()->trait_values)["ff"];

    if (os_)
        *os_ << "ff_values: " << *ff_values << endl;

    unit_assert(ff_values->size() == 20);

    for (int i=0; i<5; ++i) unit_assert(ff_values->at(i) == 1);
    for (int i=5; i<10; ++i) unit_assert(ff_values->at(i) == 0);
    for (int i=10; i<20; ++i) unit_assert(ff_values->at(i) == 0);
}


void test()
{
    test_FitnessFunction_Optimum_Polynomial();
    test_FitnessFunction_Optimum_Gaussian();
    test_FitnessFunction_TruncationSelection();
    test_FitnessFunction_TruncationSelection_2();
    test_FitnessFunction_TruncationSelection_3();
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


