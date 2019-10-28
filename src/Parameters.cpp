//
// Parameters.cpp
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
#include "boost/algorithm/string.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>


using namespace std;


void Parameters::parse(const string& arg)
{
    bool flag = false;

    size_t index_equals = arg.find('=');
    if (index_equals == string::npos)
    {
        index_equals = arg.size();
        flag = true;
    }

    string name = boost::trim_copy(arg.substr(0, index_equals));
    if (name.empty()) return;

    size_t index_value = (index_equals < arg.size()) ? index_equals + 1 : index_equals;
    size_t index_comment = arg.find('#', index_value);
    if (index_comment == string::npos) index_comment = arg.size();
    string value = flag ? "1" : boost::trim_copy(arg.substr(index_value, index_comment - index_value));

    //cout << "(name,value): (" << name << "," << value << ")\n";

    insert(make_pair(name, value));
}


void Parameters::parse_file(const string& filename)
{
    ifstream is(filename.c_str());
    if (!is) throw runtime_error(("[Parameters::parse_file()] Unable to open file " + filename).c_str());

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (buffer.empty() || buffer[0] == '#') continue;
        parse(buffer);
    }
}


ostream& operator<<(ostream& os, const Parameters& parameters)
{
    for (Parameters::const_iterator it=parameters.begin(); it!=parameters.end(); ++it)
        os << "    " << it->first << " = " << it->second << endl;

    return os;
}


