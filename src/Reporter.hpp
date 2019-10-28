//
// Reporter.hpp
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


#ifndef _REPORTER_HPP_
#define _REPORTER_HPP_


#include "Population.hpp"
#include "PopulationData.hpp"
#include "Configurable.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <ctime>


namespace bfs = boost::filesystem;


//
// Reporter 
//

///
/// interface class
///

class Reporter : public virtual Configurable
{
    public:

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation) {}

    virtual Loci loci(size_t generation_index, 
                      bool is_final_generation) const {return Loci();}

    // Configurable interface default implementation

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void initialize(const SimulatorConfig& config);
    virtual void configure(const Parameters& parameters, const Registry& registry);

    virtual ~Reporter(){}

    protected:

    bfs::path output_directory_;

    // Note: the compiler needs this because Configurable has no default
    // constructor (by design).  However, this call to Configurable(id) never
    // actually happens.  Configurable is a virtual base, so derived classes
    // must call Configurable(id) directly in their constructors, and this
    // call is ignored by compiler.
    Reporter() : Configurable("dummy_id_reporter") {}
};


typedef shared_ptr<Reporter> ReporterPtr;
typedef std::vector<ReporterPtr> ReporterPtrs;


#endif // _REPORTER_HPP_

