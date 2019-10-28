//
// RandomTest.cpp
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


#include "Random.hpp"
#include "unit.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstring>
#include <ctime>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void demo()
{
    Random::seed(std::time(0));

    if (os_) *os_ << "Random::uniform_01():\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << Random::uniform_01() << endl;
    if (os_) *os_ << endl;

    if (os_) *os_ << "Random::uniform_integer(1,10):\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << Random::uniform_integer(1,10) << endl;
    if (os_) *os_ << endl;

    if (os_) *os_ << "Random::uniform_real(5,7):\n";
    for (int i=0; i<10; i++)
        if (os_) *os_ << Random::uniform_real(5,7) << endl;
    if (os_) *os_ << endl;
}


void test_seed()
{
    if (os_) *os_ << "test_seed()\n";

    Random::seed(420); 

    if (os_) *os_ << "first pass:\n";
    vector<int> v;
    for (int i=0; i<10; i++)
    {
        v.push_back(Random::uniform_integer(1,10));
        if (os_) *os_ << v.back() << endl;
    }
    if (os_) *os_ << endl;
    
    Random::seed(420);

    if (os_) *os_ << "second pass:\n";
    vector<int> w;
    for (int i=0; i<10; i++)
    {
        w.push_back(Random::uniform_integer(1,10));
        if (os_) *os_ << w.back() << endl;
    }
    if (os_) *os_ << endl;

    for (int i=0; i<10; i++)
        unit_assert(v[i] == w[i]);
}


void test_constant_distribution()
{
    if (os_) *os_ << "test_constant_distribution()\n";

    Random::DistributionPtr d = Random::create_constant_distribution("id_dummy");

    for (int i=0; i<10; ++i)
        unit_assert(d->random_value() == 0);

    double value = 5.67;

    Parameters parameters;
    parameters.insert_name_value("value", value);

    Configurable::Registry registry;
    d->configure(parameters, registry);

    for (int i=0; i<10; ++i)
        unit_assert(d->random_value() == value);

    Parameters parameters_out = d->parameters();
    unit_assert(parameters == parameters_out);
}


void sample_mean_variance(const Random::Distribution& d, double& mean, double& variance)
{
    double sum = 0;
    double sum2 = 0;
    size_t count = 0;

    for (int i=0; i<10000; ++i)
    {
        ++count;
        double value = d.random_value();
        sum += value;
        sum2 += value*value;
    }

    mean = sum / count;
    variance = sum2/count - mean*mean;
}


void test_uniform_real_distribution() 
{
    if (os_) *os_ << "test_uniform_real_distribution()\n";

    // test default: min == 0, max == 1

    Random::DistributionPtr d = Random::create_uniform_real_distribution("id_uniform_real");

    if (os_)
    {
        *os_ << "random values: ";
        for (int i=0; i<10; ++i)
            *os_ << d->random_value() << " ";
        *os_ << endl;
    }

    double sample_mean = 0;
    double sample_variance = 0;
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    double epsilon = .01;
    unit_assert_equal(sample_mean, 0.5, epsilon);
    unit_assert_equal(sample_variance, 1./12, epsilon);

    // test configured

    double min = 1;
    double max = 13;

    Parameters parameters;
    parameters.insert_name_value("min", min);
    parameters.insert_name_value("max", max);

    Configurable::Registry registry;
    d->configure(parameters, registry);
    
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    epsilon = .1;
    unit_assert_equal(sample_mean, 7, epsilon);
    unit_assert_equal(sample_variance, 12, epsilon);

    Parameters parameters_out = d->parameters();
    unit_assert(parameters == parameters_out);
}


void test_normal_distribution() 
{
    if (os_) *os_ << "test_normal_distribution()\n";

    Random::DistributionPtr d = Random::create_normal_distribution("id_normal");

    // test default: mean == 0, variance == 1

    double sample_mean = 0;
    double sample_variance = 0;
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    double epsilon = .03;
    unit_assert_equal(sample_mean, 0.0, epsilon);
    unit_assert_equal(sample_variance, 1.0, epsilon);

    // test configured

    double mean = 9.87;
    double variance = 1.23;

    Parameters parameters;
    parameters.insert_name_value("mean", mean);
    parameters.insert_name_value("variance", variance);


    Configurable::Registry registry;
    d->configure(parameters, registry);
    
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    unit_assert_equal(sample_mean, mean, epsilon);
    unit_assert_equal(sample_variance, variance, epsilon);

    Parameters parameters_out = d->parameters();
    unit_assert(parameters == parameters_out);
}


void test_exponential_distribution() 
{
    if (os_) *os_ << "test_exponential_distribution()\n";

    Random::DistributionPtr d = Random::create_exponential_distribution("id_exponential");

    // test default: rate = 1

    double sample_mean = 0;
    double sample_variance = 0;
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    double epsilon = .08;
    unit_assert_equal(sample_mean, 1.0, epsilon);
    unit_assert_equal(sample_variance, 1.0, epsilon);

    // test configured

    double rate = 9.87;

    Parameters parameters;
    parameters.insert_name_value("rate", rate);

    Configurable::Registry registry;
    d->configure(parameters, registry);
    
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    epsilon = .002;
    unit_assert_equal(sample_mean, 1/rate, epsilon);
    unit_assert_equal(sample_variance, 1/rate/rate, epsilon);

    Parameters parameters_out = d->parameters();
    unit_assert(parameters == parameters_out);
}


void test_poisson_distribution() 
{
    if (os_) *os_ << "test_poisson_distribution()\n";

    Random::DistributionPtr d = Random::create_poisson_distribution("id_poisson");

    // test default: rate = 1

    double sample_mean = 0;
    double sample_variance = 0;
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    double epsilon = .06;
    unit_assert_equal(sample_mean, 1.0, epsilon);
    unit_assert_equal(sample_variance, 1.0, epsilon);

    // test configured

    double rate = 9.87;

    Parameters parameters;
    //parameters.insert_name_value("mean", rate);
    parameters.insert_name_value("rate", rate);

    Configurable::Registry registry;
    d->configure(parameters, registry);
    
    sample_mean_variance(*d, sample_mean, sample_variance);

    if (os_)
    {
        *os_ << "sample_mean: " << sample_mean << endl;
        *os_ << "sample_variance: " << sample_variance << endl;
    }

    epsilon = .05;
    unit_assert_equal(sample_mean, rate, epsilon);
    unit_assert_equal(sample_variance, rate, epsilon);

    Parameters parameters_out = d->parameters();
    unit_assert(parameters == parameters_out);
}


void demo_random_indices_without_replacement()
{
    if (os_) *os_ << "demo_random_indices_without_replacement()\n";

    const size_t population_size = 100;
    const size_t sample_size = 10;

    for (int i=0; i<10; ++i)
    {
        vector<size_t> indices = Random::random_indices_without_replacement(population_size, sample_size);
        if (os_)
        {
            copy(indices.begin(), indices.end(), ostream_iterator<size_t>(*os_, " "));
            *os_ << endl;
        }
    }
}


void test_discrete_distribution() 
{
    if (os_) *os_ << "test_discrete_distribution()\n";

    // test constructed

    vector<double> frequencies;
    frequencies.push_back(8);
    frequencies.push_back(4);
    frequencies.push_back(2);
    frequencies.push_back(2);

    vector<double> values;
    values.push_back(2);
    values.push_back(4);
    values.push_back(6);
    values.push_back(8);

    Random::DistributionPtr d = Random::create_discrete_distribution("id_discrete", 
            frequencies, values);

    vector<size_t> counts(10);

    size_t trial_count = 10000;
    for (size_t i=0; i<trial_count; ++i)
        ++counts[static_cast<size_t>(d->random_value())];

    vector<double> counts_normalized(10);
    transform(counts.begin(), counts.end(), counts_normalized.begin(),
              bind2nd(divides<double>(), trial_count));

    if (os_)
    {
        *os_ << "counts_normalized: ";
        copy(counts_normalized.begin(), counts_normalized.end(), 
             ostream_iterator<double>(*os_, " "));
        *os_ << endl;
    }

    const double epsilon = .03;

    unit_assert(counts_normalized[0] == 0);
    unit_assert(counts_normalized[1] == 0);
    unit_assert_equal(counts_normalized[2], .5, epsilon);
    unit_assert(counts_normalized[3] == 0);
    unit_assert_equal(counts_normalized[4], .25, epsilon);
    unit_assert(counts_normalized[5] == 0);
    unit_assert_equal(counts_normalized[6], .125, epsilon);
    unit_assert(counts_normalized[7] == 0);
    unit_assert_equal(counts_normalized[8], .125, epsilon);
    unit_assert(counts_normalized[9] == 0);

    // test configured

    Parameters parameters;
    parameters.insert_name_value("frequencies", "1.2 0.6 0.3 0.3 ");
    parameters.insert_name_value("values", "1 3 5 7 ");

    Configurable::Registry registry;
    d->configure(parameters, registry);

    Parameters parameters_out = d->parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    counts = vector<size_t>(10); // reset
    for (size_t i=0; i<trial_count; ++i)
        ++counts[static_cast<size_t>(d->random_value())];

    transform(counts.begin(), counts.end(), counts_normalized.begin(),
              bind2nd(divides<double>(), trial_count));

    if (os_)
    {
        *os_ << "counts_normalized: ";
        copy(counts_normalized.begin(), counts_normalized.end(), 
             ostream_iterator<double>(*os_, " "));
        *os_ << endl;
    }

    unit_assert(counts_normalized[0] == 0);
    unit_assert_equal(counts_normalized[1], .5, epsilon);
    unit_assert(counts_normalized[2] == 0);
    unit_assert_equal(counts_normalized[3], .25, epsilon);
    unit_assert(counts_normalized[4] == 0);
    unit_assert_equal(counts_normalized[5], .125, epsilon);
    unit_assert(counts_normalized[6] == 0);
    unit_assert_equal(counts_normalized[7], .125, epsilon);
    unit_assert(counts_normalized[8] == 0);
    unit_assert(counts_normalized[9] == 0);
}


void test_neutral_frequency_distribution() 
{
    if (os_) *os_ << "test_neutral_frequency_distribution()\n";

    Random::DistributionPtr d = Random::create_neutral_frequency_distribution("id_neutral_frequency");

    const unsigned int sample_size = 100;

    Parameters parameters;
    parameters.insert_name_value("sample_size", sample_size);

    Configurable::Registry registry;
    d->configure(parameters, registry);

    Parameters parameters_out = d->parameters();

    if (os_)
    {
        *os_ << "parameters:\n" << parameters << endl;
        *os_ << "parameters_out:\n" << parameters_out << endl;
    }

    unit_assert(parameters == parameters_out);

    if (os_)
    {
        *os_ << "samples: ";
        for (int i=0; i<10; i++)
            *os_ << d->random_value() << " ";
        *os_ << endl;
    }

    const size_t trial_count = 10000;
    vector<size_t> counts(sample_size);
    for (size_t i=0; i<trial_count; ++i)
        ++counts[static_cast<size_t>(d->random_value() * sample_size)];

    vector<double> counts_normalized(sample_size);
    transform(counts.begin(), counts.end(), counts_normalized.begin(),
              bind2nd(divides<double>(), trial_count));

    vector<double> ratios(sample_size);
    transform(counts.begin(), counts.end(), ratios.begin(),
              bind2nd(divides<double>(), counts[1]));

    if (os_)
    {
        *os_ << "trial count: " << trial_count << endl;
        *os_ << "sample_size: " << sample_size << endl;

        *os_ << "counts_normalized: ";
        copy(counts_normalized.begin(), counts_normalized.end(), 
             ostream_iterator<double>(*os_, " "));
        *os_ << endl;

        *os_ << "ratios: ";
        copy(ratios.begin(), ratios.end(), 
             ostream_iterator<double>(*os_, " "));
        *os_ << endl;
    }


    const double epsilon = .05;

    unit_assert(counts_normalized[0] == 0);
    unit_assert_equal(ratios[2], 1./2, epsilon);
    unit_assert_equal(ratios[3], 1./3, epsilon);
    unit_assert_equal(ratios[4], 1./4, epsilon);
    unit_assert_equal(ratios[5], 1./5, epsilon);
}


void test()
{
    demo();
    test_seed();
    test_constant_distribution();
    test_uniform_real_distribution();
    test_normal_distribution();
    test_exponential_distribution();
    test_poisson_distribution();
    demo_random_indices_without_replacement();
    test_discrete_distribution();
    test_neutral_frequency_distribution();
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


