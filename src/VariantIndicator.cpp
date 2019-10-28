//
// VariantIndicator.cpp
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


#include "VariantIndicator.hpp"


using namespace std;


//
// VariantIndicator
//


void VariantIndicator::write_file(const std::string& filename) const
{
    ofstream os(filename.c_str());
    os << class_name() << " " << object_id() << endl << parameters() << endl;
}


unsigned int VariantIndicator::mutate(unsigned int old_chunk_id, const Locus& locus, unsigned int value)
{
    throw runtime_error("[VariantIndicator::mutate()] Not implemented.");
}


string VariantIndicator::class_name() const
{
    cerr << "[VariantIndicator] Warning: virtual class_name() has not been defined in derived class.\n";
    return "VariantIndicator";
}


Parameters VariantIndicator::parameters() const
{
    cerr << "[VariantIndicator] Warning: virtual parameters() has not been defined in derived class.\n";
    return Parameters();
}


void VariantIndicator::configure(const Parameters& parameters, const Registry& registry)
{
    cerr << "[VariantIndicator] Warning: virtual configure() has not been defined in derived class.\n";
}


