//
// Parameters.hpp
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


#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_


#include "boost/lexical_cast.hpp"
#include <string>
#include <map>
#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <iterator>
#include <stdexcept>


class Parameters : public std::multimap<std::string, std::string>
{
    public:

    void parse(const std::string& arg);
    void parse_file(const std::string& filename);

    // insertion with lexical cast to string

    template <typename result_type>
    void insert_name_value(const std::string& name, const result_type& value)
    {
        insert(make_pair(name, boost::lexical_cast<std::string>(value)));
    }

    template <typename result_type>
    void insert_name_value_vector(const std::string& name, const std::vector<result_type>& values)
    {
        std::ostringstream oss;
        copy(values.begin(), values.end(), std::ostream_iterator<result_type>(oss, " "));
        insert(make_pair(name, oss.str()));
    }

    //
    // value(): retrieve (first) value with type conversion
    //
    // overloaded versions have different semantics:
    //  1) value(name)          : if name not found, throw runtime_error
    //  2) value(name, default) : if name not found, return default
    //

    template <typename result_type>
    result_type value(const std::string& name) const
    {
        accessed.insert(name);
        if (count(name) == 0) 
            throw std::runtime_error("[Parameters::value()] Required parameter not found: " + name);
        std::string value = lower_bound(name)->second; // return first value
        return static_cast<result_type>(boost::lexical_cast<double>(value));
    }

    template <typename result_type>
    result_type value(const std::string& name, result_type default_value) const
    {
        accessed.insert(name);
        if (count(name) == 0) return default_value;
        std::string value = lower_bound(name)->second; // return first value
        return static_cast<result_type>(boost::lexical_cast<double>(value));
    }

    //
    // values(): retrieve vector of all (multi-mapped) values, with type conversion (may be empty)
    //           

    template <typename result_type>
    std::vector<result_type> values(const std::string& name) const
    {
        accessed.insert(name);
        std::vector<result_type> result;
        std::pair<const_iterator, const_iterator> range = equal_range(name);
        for (Parameters::const_iterator it=range.first; it!=range.second; ++it)
            result.push_back(static_cast<result_type>(boost::lexical_cast<double>(it->second)));
        return result;
    }

    //
    // value_vector(): parses/separates each value via istream iteration, with type conversion
    //
    // overloaded versions have different semantics:
    //  1) value_vector(name)          : if name not found, throw runtime_error
    //  2) value_vector(name, default) : if name not found, return default
    //

    template <typename result_type>
    std::vector<result_type> value_vector_aux(const std::string& name) const
    {
        accessed.insert(name);
        std::vector<result_type> result;
        std::pair<const_iterator, const_iterator> range = equal_range(name);
        for (Parameters::const_iterator it=range.first; it!=range.second; ++it)
        {
            std::istringstream iss(it->second);

            std::istream_iterator<double> jt = std::istream_iterator<double>(iss);
            for (; jt!=std::istream_iterator<double>(); ++jt)
                result.push_back(static_cast<result_type>(*jt));
        }
        return result;
    }

    template <typename result_type>
    std::vector<result_type> value_vector(const std::string& name) const
    {
        if (count(name) == 0) 
            throw std::runtime_error("[Parameters::value_vector()] Required parameter not found: " + name);
        return value_vector_aux<result_type>(name);
    }

    template <typename result_type>
    std::vector<result_type> value_vector(const std::string& name, const std::vector<result_type>& default_value) const
    {
        if (count(name) == 0) return default_value;
        return value_vector_aux<result_type>(name);
    }

    //
    // track parameters that have been retrived via value*()
    //

    mutable std::set<std::string> accessed;
};


//
// template specializations for string
//


template <>
inline void Parameters::insert_name_value<std::string>(const std::string& name, const std::string& value)
{
    insert(make_pair(name, value));
}


template <> 
inline std::string Parameters::value<std::string>(const std::string& name) const
{
    accessed.insert(name);
    if (count(name) == 0) 
        throw std::runtime_error(("[Parameters::value()] Required parameter not found: " + name).c_str());
    return lower_bound(name)->second;
}


template <> 
inline std::string Parameters::value<std::string>(const std::string& name, std::string default_value) const
{
    accessed.insert(name);
    if (count(name) == 0) return default_value;
    return lower_bound(name)->second;
}


template <>
inline std::vector<std::string> Parameters::values<std::string>(const std::string& name) const
{
    accessed.insert(name);
    std::vector<std::string> result;
    std::pair<const_iterator, const_iterator> range = equal_range(name);
    for (Parameters::const_iterator it=range.first; it!=range.second; ++it)
        result.push_back(it->second);
    return result;
}


template <>
inline std::vector<std::string> Parameters::value_vector_aux<std::string>(const std::string& name) const
{
    accessed.insert(name);
    std::vector<std::string> result;
    std::pair<const_iterator, const_iterator> range = equal_range(name);
    for (Parameters::const_iterator it=range.first; it!=range.second; ++it)
    {
        std::istringstream iss(it->second);
        copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
             back_inserter(result));
    }
    return result;
}


//
// notes: 
//
// 1) to allow integers in scientific notation, we do lexical_cast<double> first, 
//    followed by static_cast<numeric_type>
// 
// 2) explicit template specializations of member functions must be declared in the namespace
//    of the enclosing class (i.e. outside the class)
//
// 3) the specialization acts like a normal function: either define as inline in header, or 
//    declare in header and define in cpp (else linker errors)
//


std::ostream& operator<<(std::ostream& os, const Parameters& parameters);


#endif // _PARAMETERS_HPP_

