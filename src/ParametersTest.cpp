//
// ParametersTest.cpp
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


#include "Parameters.hpp"
#include "unit.hpp"
#include "boost/math/constants/constants.hpp"
#include <iostream>
#include <cstring>
#include <iterator>
#include <set>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_parse()
{
    if (os_) *os_ << "test_parse()\n";

    Parameters p;
    p.parse("  Simba=red");
    p.parse("   Gabi = black  ");
    p.parse("\t \r Gizmo \t\t = \n  black and white   ");
    p.parse("  \n   Lucky  =\tblack  \n\n\r");
    p.parse("  verbose   ");
    p.parse("pi=3.14159");
    p.parse("commented_value =  123 # this stuff should be ignored");
    p.parse(" population_size  =  1.23e6"); // scientific notation for integers

    if (os_) *os_ << p;

    unit_assert(p.value<string>("Simba") == "red");
    unit_assert(p.value<string>("Gabi") == "black");
    unit_assert(p.value<string>("Gizmo") == "black and white");
    unit_assert(p.value<string>("Lucky") == "black");
    unit_assert(p.value<int>("verbose") == 1);
    unit_assert(p.value<bool>("verbose") == true);
    unit_assert(p.value<double>("pi") == 3.14159);
    unit_assert(p.value<double>("e", 2.71828) == 2.71828); // specified default value
    unit_assert(p.value<string>("e", "f") == "f"); // specified default value
    unit_assert(p.value<int>("commented_value") == 123);
    unit_assert(p.value<unsigned int>("population_size") == 1230000);

    bool caught = false;
    try
    {
        double temp = p.value<double>("e");
        (void)temp;
    }
    catch (runtime_error& e)
    {
        if (os_) *os_ << "Caught: " << e.what() << "\n\n";;
        caught = true;
    }
    unit_assert(caught);
}


void test_insert()
{
    if (os_) *os_ << "test_insert()\n";

    Parameters p;
    p.insert(make_pair("one", "1"));
    p.insert_name_value("two", 2);
    p.insert_name_value("three", string("3"));

    if (os_) *os_ << p << endl;

    unit_assert(p.value<int>("one") == 1);
    unit_assert(p.value<int>("two") == 2);
    unit_assert(p.value<int>("three") == 3);
}


void test_values()
{
    if (os_) *os_ << "test_values()\n";

    const double pi = boost::math::constants::pi<double>();
    const double e = boost::math::constants::e<double>();

    Parameters p;
    p.insert_name_value("favorite_numbers", 10);
    p.insert_name_value("favorite_numbers", 27);
    p.insert_name_value("favorite_numbers", 72);
    p.insert_name_value("favorite_transcendentals", pi);
    p.insert_name_value("favorite_transcendentals", e);

    vector<int> favorite_numbers = p.values<int>("favorite_numbers");
    vector<double> favorite_transcendentals = p.values<double>("favorite_transcendentals");
    vector<int> favorite_transcendentals_as_ints = p.values<int>("favorite_transcendentals");

    if (os_)
    {
        *os_ << "favorite_numbers: ";
        copy(favorite_numbers.begin(), favorite_numbers.end(), ostream_iterator<int>(*os_, " "));
        *os_ << endl;

        *os_ << "favorite_transcendentals: ";
        copy(favorite_transcendentals.begin(), favorite_transcendentals.end(), ostream_iterator<double>(*os_, " "));
        *os_ << endl;

        *os_ << "favorite_transcendentals_as_ints: ";
        copy(favorite_transcendentals_as_ints.begin(), favorite_transcendentals_as_ints.end(), ostream_iterator<double>(*os_, " "));
        *os_ << endl;
    }

    unit_assert(favorite_numbers.size() == 3);
    unit_assert(favorite_transcendentals.size() == 2);
    unit_assert(favorite_transcendentals_as_ints.size() == 2);

    // note: order of values not necessarily preserved in multimap

    set<int> favorite_numbers_unordered;
    copy(favorite_numbers.begin(), favorite_numbers.end(), inserter(favorite_numbers_unordered, favorite_numbers_unordered.begin()));
    unit_assert(favorite_numbers_unordered.count(10));
    unit_assert(favorite_numbers_unordered.count(27));
    unit_assert(favorite_numbers_unordered.count(72));
    unit_assert(!favorite_numbers_unordered.count(0));

    const double epsilon = 1e-5;

    if (favorite_transcendentals[0] > 3)
    {
        unit_assert_equal(favorite_transcendentals[0], pi, epsilon);
        unit_assert_equal(favorite_transcendentals[1], e, epsilon);
    }
    else
    {
        unit_assert_equal(favorite_transcendentals[1], pi, epsilon);
        unit_assert_equal(favorite_transcendentals[0], e, epsilon);
    }

    set<int> favorite_transcendentals_as_ints_unordered;
    copy(favorite_transcendentals_as_ints.begin(), favorite_transcendentals_as_ints.end(), inserter(favorite_transcendentals_as_ints_unordered, favorite_transcendentals_as_ints_unordered.begin()));
    unit_assert(favorite_transcendentals_as_ints_unordered.count(2));
    unit_assert(favorite_transcendentals_as_ints_unordered.count(3));
}


void test_value_vector()
{
    if (os_) *os_ << "test_value_vector()\n";

    vector<string> cat_names;
    cat_names.push_back("Gabi");
    cat_names.push_back("Simba");
    cat_names.push_back("Gizmo");
    cat_names.push_back("Lucky");

    vector<string> more_cat_names;
    more_cat_names.push_back("Felix");
    more_cat_names.push_back("Garfield");
    more_cat_names.push_back("Heathcliff");

    Parameters parameters;
    parameters.insert_name_value("even", "0 2 4");
    parameters.insert_name_value("even", "6 8");
    parameters.insert_name_value("transcendentals", "2.71828 3.14159");
    parameters.insert_name_value_vector("names", cat_names);
    parameters.insert_name_value_vector("names", more_cat_names);

    vector<int> default_odd;
    for (int i=1; i<10; i+=2) default_odd.push_back(i);

    vector<int> even_numbers = parameters.value_vector<int>("even");
    vector<int> even_numbers_2 = parameters.value_vector<int>("even", vector<int>());
    vector<int> odd_numbers = parameters.value_vector<int>("odd", default_odd);
    vector<double> transcendentals = parameters.value_vector<double>("transcendentals");
    vector<unsigned int> transcendentals_as_ints = parameters.value_vector<unsigned int>("transcendentals");
    vector<string> names = parameters.value_vector<string>("names");

    if (os_)
    {
        *os_ << "even_numbers: ";
        copy(even_numbers.begin(), even_numbers.end(), ostream_iterator<int>(*os_, " "));
        *os_ << endl;

        *os_ << "even_numbers_2: ";
        copy(even_numbers_2.begin(), even_numbers_2.end(), ostream_iterator<int>(*os_, " "));
        *os_ << endl;

        *os_ << "odd_numbers: ";
        copy(odd_numbers.begin(), odd_numbers.end(), ostream_iterator<int>(*os_, " "));
        *os_ << endl;

        *os_ << "transcendentals: ";
        copy(transcendentals.begin(), transcendentals.end(), ostream_iterator<double>(*os_, " "));
        *os_ << endl;

        *os_ << "transcendentals_as_ints: ";
        copy(transcendentals_as_ints.begin(), transcendentals_as_ints.end(), ostream_iterator<unsigned int>(*os_, " "));
        *os_ << endl;

        *os_ << "names: ";
        copy(names.begin(), names.end(), ostream_iterator<string>(*os_, " "));
        *os_ << endl;
    }

    unit_assert(even_numbers.size() == 5);
    unit_assert(even_numbers[0] == 0);
    unit_assert(even_numbers[4] == 8);
    unit_assert(even_numbers_2.size() == 5);
    unit_assert(even_numbers_2[0] == 0);
    unit_assert(even_numbers_2[4] == 8);
    unit_assert(odd_numbers.size() == 5);
    unit_assert(transcendentals.size() == 2);
    unit_assert(transcendentals_as_ints.size() == 2);
    unit_assert(transcendentals_as_ints[0] == 2);
    unit_assert(transcendentals_as_ints[1] == 3);
    unit_assert(names.size() == 7);
    unit_assert(names[0] == "Gabi");
    unit_assert(names[1] == "Simba");
    unit_assert(names[2] == "Gizmo");
    unit_assert(names[3] == "Lucky");
    unit_assert(names[4] == "Felix");
    unit_assert(names[5] == "Garfield");
    unit_assert(names[6] == "Heathcliff");
}


void test_accessed()
{
    if (os_) *os_ << "test_accessed()\n";

    Parameters parameters;
    parameters.insert_name_value("cat", "gabi");
    parameters.insert_name_value("dog", "fido");
    parameters.insert_name_value("horse", "ed");
    parameters.insert_name_value("e", 2.71828);
    parameters.insert_name_value("pi", 3);

    parameters.value<string>("cat");
    parameters.value<int>("pi");
    parameters.value<double>("e");

    if (os_)
    {
        *os_ << "accessed: ";
        copy(parameters.accessed.begin(), parameters.accessed.end(), ostream_iterator<string>(*os_, " "));
        *os_ << endl;
    }

    unit_assert(parameters.accessed.size() == 3);
    unit_assert(parameters.accessed.count("cat"));
    unit_assert(parameters.accessed.count("pi"));
    unit_assert(parameters.accessed.count("e"));
}


void test()
{
    test_parse();
    test_insert();
    test_values();
    test_value_vector();
    test_accessed();
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


