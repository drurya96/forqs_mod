//
// DataVector.cpp
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
#include <iostream>
#include <iterator>
#include <numeric>
#include <fstream>
#include <stdexcept>


using namespace std;


double DataVector::sum() const
{
    return accumulate(begin(), end(), 0.0);
}



double DataVector::mean() const
{
    if (empty()) return 0;
    return accumulate(begin(), end(), 0.0)/size();
}


double DataVector::variance() const
{
    double m = mean();
    DataVector d2 = *this;
    d2 *= *this;
    return d2.mean() - m*m;
}


DataVectorPtr DataVector::cdf() const
{
    DataVectorPtr result(new DataVector(this->size()));
    partial_sum(this->begin(), this->end(), result->begin());    
    return result;
}


bool DataVector::all_zero() const
{
    for (const_iterator it=begin(); it!=end(); ++it)
        if (*it != 0) return false;

    return true;
}


void DataVector::write_file(const std::string& filename) const
{
    ofstream os(filename.c_str());
    if (!os)
        throw runtime_error(("[DataVector::write_file()] Unable to open file " + filename).c_str());

    copy(begin(), end(), ostream_iterator<double>(os, "\n"));
}


void DataVector::operator*=(const DataVector& that)
{
    if (that.size() != this->size())
        throw runtime_error("[DataVector::operator*=] Vector size mismatch.");

    const_iterator jt = that.begin();
    for (iterator it=begin(); it!=end(); ++it, ++jt)
        *it *= *jt;
}


ostream& operator<<(ostream& os, const DataVector& v)
{
    copy(v.begin(), v.end(), ostream_iterator<double>(os, " "));
    return os;
}


