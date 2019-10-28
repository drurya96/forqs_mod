//
// QuantitativeTraitTest.cpp
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


#include "QuantitativeTrait.hpp"
#include "unit.hpp"
#include <iostream>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_map_get()
{
    DataVectorPtr data(new DataVector);
    string id_good("id_good");
    string id_bad("id_bad");
    TraitValueMap tvm;
    tvm[id_good] = data;

    DataVectorPtr retrieved = tvm.get(id_good); // ok
    unit_assert(retrieved.get() == data.get());

    // try to retrieve data that isn't in the map

    bool caught = false;
    try
    {
        DataVectorPtr dummy = tvm.get(id_bad);
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught exception:\n" << e.what() << endl;
        caught = true;
    }
    unit_assert(caught);

    // try to retrieve null data vector

    DataVectorPtr null;
    tvm[id_bad] = null;
    caught = false;
    try
    {
        DataVectorPtr dummy = tvm.get(id_bad);
    }
    catch (exception& e)
    {
        if (os_) *os_ << "caught exception:\n" << e.what() << endl;
        caught = true;
    }
    unit_assert(caught);
}


void test()
{
    test_map_get();
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

