//
// Random.hpp
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


#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_


#include "Configurable.hpp"
#include "shared_ptr.hpp"


//
// simple wrapper for Boost.Random
//


class Random
{
    public:
    
    // set seed
    static void seed(unsigned int value);

    // return random integer N with a <= N <= b
    static int uniform_integer(int a, int b);

    // return random long integer N with a <= N <= b
    static long uniform_long(long a, long b);

    // return random double in [a,b)
    static double uniform_real(double a, double b);

    // return random double in [0,1)
    static double uniform_01();

    // return 1 with probability p, else 0
    static int bernoulli(double p = .5);

    // return random sample of indices without replacement
    static std::vector<size_t> random_indices_without_replacement(size_t population_size, size_t sample_size);

    // distributions

    class Distribution;
    typedef shared_ptr<Distribution> DistributionPtr;

    static DistributionPtr create_constant_distribution(const std::string& id, 
                                                        double value = 0);

    static DistributionPtr create_uniform_real_distribution(const std::string& id, 
                                                            double min = 0, double max = 1);

    static DistributionPtr create_normal_distribution(const std::string& id, 
                                                      double mean = 0, double variance = 1);

    static DistributionPtr create_exponential_distribution(const std::string& id, 
                                                           double rate = 1);

    static DistributionPtr create_poisson_distribution(const std::string& id, 
                                                       double rate = 1);

    static DistributionPtr create_discrete_distribution(const std::string& id, 
                                                        std::vector<double> frequencies = std::vector<double>(),
                                                        std::vector<double> values = std::vector<double>());

    static DistributionPtr create_neutral_frequency_distribution(const std::string& id, 
                                                                 unsigned int sample_size = 0);
};


///
/// \defgroup Distributions Random Distributions
///
/// Configurable probability distributions
///


//
// Random::Distribution
//


///
/// interface class
///


class Random::Distribution : public Configurable
{
    public:

    Distribution(const std::string& id) : Configurable(id) {}

    virtual double random_value() const = 0;

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
};


#endif // _RANDOM_HPP_

