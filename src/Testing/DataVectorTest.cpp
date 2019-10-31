//
// DataVectorTest.cpp
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


#include "DataVector.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_cdf()
{
    if (os_) *os_ << "test_cdf()\n";

    DataVectorPtr a(new DataVector(5, 1.0));

    if (os_) *os_ << *a << endl;

    DataVectorPtr b = a->cdf();

    if (os_) *os_ << *b << endl;

    for (size_t i=0; i<5; ++i)
        unit_assert(b->at(i) == i+1);
}


void test_mean()
{
    DataVector d;
    for (int i=1; i<=9; ++i)
        d.push_back(i);

    unit_assert(d.sum() == 45.0);
    unit_assert(d.mean() == 5.0);
}


void test_multiplication_update()
{
    DataVector d;
    for (int i=0; i<10; ++i) d.push_back(i);

    d *= d;

    for (int i=0; i<10; ++i) unit_assert(d[i] == i*i);
}


void test_variance()
{
    DataVector d;
    for (size_t i=0; i<7; ++i) d.push_back(1);
    for (size_t i=0; i<3; ++i) d.push_back(0);
    const double epsilon = 1e-6;
    unit_assert_equal(d.variance(), .7 * .3, epsilon);
}


void test_all_zero()
{
    DataVector d;
    unit_assert(d.all_zero());
    for (size_t i=0; i<3; ++i) d.push_back(0);
    unit_assert(d.all_zero());
    d.push_back(1);
    unit_assert(!d.all_zero());
}


void test()
{
    test_cdf();
    test_mean();
    test_multiplication_update();
    test_variance();
    test_all_zero();
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


