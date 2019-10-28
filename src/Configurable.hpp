//
// Configurable.hpp
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


#ifndef _CONFIGURABLE_HPP_
#define _CONFIGURABLE_HPP_


#include "Parameters.hpp"
#include "shared_ptr.hpp"
#include <string>
#include <set>
#include <stdexcept>


struct SimulatorConfig;


class Configurable
{
    public:

    class Registry;

    Configurable(const std::string& id) : id_(id) {}
    const std::string& object_id() const {return id_;}

    virtual std::string class_name() const = 0;
    virtual Parameters parameters() const = 0;
    virtual void configure(const Parameters& parameters, const Registry& registry) = 0;
    virtual void initialize(const SimulatorConfig& config) {}

    void write_configuration(std::ostream& os, std::set<std::string>& ids_written) const;
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const {}

    virtual ~Configurable() {}

    private:

    std::string id_;
};


typedef shared_ptr<Configurable> ConfigurablePtr;
typedef std::vector<ConfigurablePtr> ConfigurablePtrs;


struct Locus;


class Configurable::Registry : public std::map<std::string, ConfigurablePtr> 
{
    public:

    template <typename result_type>
    shared_ptr<result_type> get(const std::string& name) const
    {
        if (!count(name))
            throw std::runtime_error(("[Configurable::Registry] Object id \"" + name + "\" requested, but not found in registry.").c_str());

        shared_ptr<result_type> result = dynamic_pointer_cast<result_type>(at(name));
        if (!result.get()) // dynamic_pointer_cast doesn't throw
            throw std::runtime_error(("[Configurable::Registry] Unable to convert object " + name).c_str());

        return result;
    }

    template <typename result_type>
    shared_ptr<result_type> get(const std::string& name, std::nothrow_t) const
    {
        if (!count(name)) return shared_ptr<result_type>();
        return dynamic_pointer_cast<result_type>(at(name));
    }
};


template<> shared_ptr<Locus> Configurable::Registry::get(const std::string& name) const;
template<> shared_ptr<Locus> Configurable::Registry::get(const std::string& name, std::nothrow_t) const;


#endif // _CONFIGURABLE_HPP_


