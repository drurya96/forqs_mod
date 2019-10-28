//
// unit.hpp
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


#ifndef _UNIT_HPP_
#define _UNIT_HPP_


#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>


//
// These are assertion macros for unit testing.  They throw a runtime_error 
// exception on failure, instead of calling abort(), allowing the application
// to recover and return an appropriate error value to the shell.
//
// unit_assert(x):                             asserts x is true
// unit_assert_equal(x, y, epsilon):           asserts x==y, within epsilon
// unit_assert_matrices_equal(A, B, epsilon):  asserts A==B, within epsilon
//


inline std::string unit_assert_message(const char* filename, int line, const char* expression)
{
    std::ostringstream oss;
    oss << "[" << filename << ":" << line << "] Assertion failed: " << expression;
    return oss.str();
}

inline std::string unit_assert_equal_message(const char* filename, int line, double x, double y, double epsilon)
{
    std::ostringstream oss;
    oss.precision(10);
    oss << "[" << filename << ":" << line << "] Assertion failed: |" << x << " - " << y << "| < " << epsilon;
    return oss.str();
}

inline std::string unit_assert_exception_message(const char* filename, int line, const char* expression, const std::string& exception)
{
    std::ostringstream oss;
    oss << "[" << filename << ":" << line << "] Assertion failed to throw \"" << exception << "\": " << expression;
    return oss.str();
}

#define unit_assert(x) \
    (!(x) ? throw std::runtime_error(unit_assert_message(__FILE__, __LINE__, #x)) : 0) 


#define unit_assert_equal(x, y, epsilon) \
    (!(fabs((x)-(y)) <= (epsilon)) ? throw std::runtime_error(unit_assert_equal_message(__FILE__, __LINE__, (x), (y), (epsilon))) : 0)


#define unit_assert_throws(x, exception) \
    { \
        bool threw = false; \
        try { (x); } \
        catch (exception&) \
        { \
            threw = true; \
        } \
        if (!threw) \
            throw std::runtime_error(unit_assert_exception_message(__FILE__, __LINE__, #x, #exception)); \
    }


#define unit_assert_throws_what(x, exception, whatStr) \
    { \
        bool threw = false; \
        try { (x); } \
        catch (exception& e) \
        { \
            if (e.what() == std::string(whatStr)) \
                threw = true; \
            else \
                throw std::runtime_error(unit_assert_exception_message(__FILE__, __LINE__, #x, std::string(#exception)+" "+(whatStr)+"\nBut a different exception was thrown: ")+(e.what())); \
        } \
        if (!threw) \
            throw std::runtime_error(unit_assert_exception_message(__FILE__, __LINE__, #x, std::string(#exception)+" "+(whatStr))); \
    }


#define unit_assert_matrices_equal(A, B, epsilon) \
    unit_assert(boost::numeric::ublas::norm_frobenius((A)-(B)) < (epsilon))


#define unit_assert_vectors_equal(A, B, epsilon) \
    unit_assert(boost::numeric::ublas::norm_2((A)-(B)) < (epsilon))


#endif // _UNIT_HPP_

