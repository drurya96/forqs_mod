//
// VariantIndicator.hpp
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


#ifndef _VARIANTINDICATOR_HPP_
#define _VARIANTINDICATOR_HPP_


#include "Configurable.hpp"
#include "Locus.hpp"
#include "Reporter.hpp"
#include "Simulator.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"


//
// VariantIndicator
//

/// 
/// interface class
///

class VariantIndicator : public virtual Configurable
{
    public:

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const = 0;
    virtual void write_file(const std::string& filename) const;
    virtual unsigned int mutate(unsigned int old_chunk_id, const Locus& locus, unsigned int value); // returns new chunk id
    virtual ~VariantIndicator() {}

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    protected:

    // Note: the compiler needs this because Configurable has no default
    // constructor (by design).  However, this call to Configurable(id) never
    // actually happens.  Configurable is a virtual base, so derived classes
    // must call Configurable(id) directly in their constructors, and this
    // call is ignored by compiler.
    VariantIndicator() : Configurable("dummy_id_variant_indicator") {}
};


typedef shared_ptr<VariantIndicator> VariantIndicatorPtr;
typedef std::vector<VariantIndicatorPtr> VariantIndicatorPtrs;


#endif // _VARIANTINDICATOR_HPP_

