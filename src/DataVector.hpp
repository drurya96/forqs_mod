//
// DataVector.hpp
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


#ifndef _DATAVECTOR_HPP_
#define _DATAVECTOR_HPP_


#include "shared_ptr.hpp"
#include <iosfwd>
#include <vector>
#include <map>
#include <stdexcept>
#include <string>


class DataVector;
typedef shared_ptr<DataVector> DataVectorPtr;
typedef std::vector<DataVectorPtr> DataVectorPtrs;


class DataVector : public std::vector<double>
{
    public:

    DataVector(size_t n = 0, double value = 0)
    :   std::vector<double>(n, value)
    {}

    // convenience functions

    double sum() const;
    double mean() const;
    double mean_nonzero() const;
    double variance() const;
    DataVectorPtr cdf() const; // note: allocates a new DataVector
    bool all_zero() const;

    void write_file(const std::string& filename) const;

    void operator*=(const DataVector& that); // element-wise multiplication update

    // future:  single precomputation of sum, sum of squares, mean, variance etc.
};


std::ostream& operator<<(std::ostream& os, const DataVector& v);


class TraitValueMap : public std::map<std::string, DataVectorPtr> // map QT id -> trait_values
{
    public:

    DataVectorPtr get(const std::string& qtid) const // returns valid pointer, or throws
    {
        if (!count(qtid))
            throw std::runtime_error(("[TraitValueMap] Quantitative trait id " + qtid + " not found.").c_str());

        DataVectorPtr result = at(qtid);

        if (!result.get())
            throw std::runtime_error(("[TraitValueMap] Null data for id " + qtid).c_str());

        return result;
    }
};


typedef shared_ptr<TraitValueMap> TraitValueMapPtr;


#endif //  _DATAVECTOR_HPP_

