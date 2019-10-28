//
// muparser_test.cpp
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


#include "muparser/muParser.h"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <stdexcept>


using namespace mu;
using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


struct Data
{
    size_t sample_size;
    double threshold;
    vector<double> trait_values;
    vector<double> male;

    Data(size_t _sample_size)
    :   sample_size(_sample_size),
        threshold(sample_size * .5),
        trait_values(sample_size),
        male(sample_size)
    {
        for (size_t i=0; i<sample_size; ++i)
        {
            trait_values[i] = i;
            male[i] = i%2;
        }
    }
};


ostream& operator<<(ostream& os, const vector<double>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<double>(os, " "));
    return os;
}


ostream& operator<<(ostream& os, const Data& data)
{
    os << "sample_size: " << data.sample_size << endl;
    os << "threshold: " << data.threshold << endl;
    os << "trait_values: " << data.trait_values << endl;;
    os << "male: " << data.male << endl;
    return os;
}


double calculate_fitness(double trait_value, double threshold, double male)
{
    return (trait_value > threshold) * male;
}


vector<double> calculate_fitness_values_1(const Data& data)
{
    // normal C++ calculation

    vector<double> fitness(data.sample_size);

    for (size_t i=0; i<data.sample_size; ++i)
    {
        fitness[i] = calculate_fitness(data.trait_values[i], data.threshold, data.male[i]);
    }

    return fitness;
}


vector<double> calculate_fitness_values_2(const Data& data)
{
    // using muparser in a loop, updating variables at each iteration

    double trait_value = 0;
    double male = 0;
    double threshold = data.threshold;

    Parser p;
    p.DefineVar("v", &trait_value);
    p.DefineVar("m", &male);
    p.DefineConst("t", threshold);
    p.SetExpr("(v>t)*m");

    vector<double> fitness(data.sample_size);

    for (size_t i=0; i<data.sample_size; ++i)
    {
        trait_value = data.trait_values[i];
        male = data.male[i];
        fitness[i] = p.Eval();
    }

    return fitness;
}


vector<double> calculate_fitness_values_3(Data& data)
{
    // using muparser array evaluation

    Parser p;
    p.DefineVar("v", &data.trait_values[0]);
    p.DefineVar("m", &data.male[0]);
    p.DefineConst("t", data.threshold);
    p.SetExpr("(v>t)*m");

    vector<double> fitness(data.sample_size);
    p.Eval(&fitness[0], fitness.size());
    return fitness;
}


void test_muparser()
{
    if (os_) *os_ << "test_muparser()\n";

    try
    {
        Data data(10);
        if (os_) *os_ << data << endl;

        vector<double> fitness1 = calculate_fitness_values_1(data);
        if (os_) *os_ << "fitness1: " << fitness1 << endl;

        vector<double> fitness2 = calculate_fitness_values_2(data);
        if (os_) *os_ << "fitness2: " << fitness2 << endl;

        vector<double> fitness3 = calculate_fitness_values_3(data);
        if (os_) *os_ << "fitness3: " << fitness3 << endl;

        for (size_t i=0; i<data.sample_size; ++i)
        {
            unit_assert(fitness1[i] == fitness2[i]);
            unit_assert(fitness1[i] == fitness3[i]);
        }
    }
    catch (Parser::exception_type &e)
    {
        if (os_) *os_ << e.GetMsg() << std::endl;
    }
}


void test()
{
    test_muparser();
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


