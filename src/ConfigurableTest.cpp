//
// ConfigurableTest.cpp
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


#include "Configurable.hpp"
#include "Locus.hpp"
#include "unit.hpp"
#include <iostream>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


class A : public Configurable 
{
    public:
    A() : Configurable("A") {}
    virtual std::string class_name() const {return "A";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}
};

typedef shared_ptr<A> APtr;


class B : public Configurable
{
    public:
    B() : Configurable("B") {}
    virtual std::string class_name() const {return "B";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}
};

typedef shared_ptr<B> BPtr;


void test_registry()
{
    if (os_) *os_ << "test_registry()\n";

    Configurable::Registry registry;

    registry["a"] = APtr(new A);
    registry["b"] = BPtr(new B);

    APtr a = registry.get<A>("a");
    unit_assert(a.get());

    BPtr b = registry.get<B>("b");
    unit_assert(b.get());

    // retrieval error

    bool caught = false;

    try
    {
        BPtr b_wrong_id = registry.get<B>("c");
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught: " << e.what() << endl;
        caught = true;
    }

    unit_assert(caught);

    // conversion error

    caught = false;

    try
    {
        BPtr b_bad_conversion = registry.get<B>("a");
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught: " << e.what() << endl;
        caught = true;
    }

    unit_assert(caught);

    // nothrow

    caught = false;

    try
    {
        BPtr b_bad_conversion = registry.get<B>("a", std::nothrow);
        unit_assert(!b_bad_conversion.get());
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught: " << e.what() << endl;
        caught = true;
    }

    unit_assert(!caught);
}


void test_get_from_locus_list()
{
    LocusListPtr locus_list(new LocusList("qtls"));
    locus_list->push_back(Locus("blah", 0, 123));
    locus_list->push_back(Locus("goo", 1, 234));

    Configurable::Registry registry;
    registry["qtls"] = locus_list;

    LocusPtr blah = registry.get<Locus>("qtls[0]");
    unit_assert(blah->chromosome_pair_index == 0 && blah->position == 123);
    LocusPtr goo = registry.get<Locus>("qtls[1]");
    unit_assert(goo->chromosome_pair_index == 1 && goo->position == 234);

    LocusPtr bad = registry.get<Locus>("qtls[2]", std::nothrow);
    unit_assert(!bad.get());
}


void test()
{
    test_registry();     
    test_get_from_locus_list();
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


