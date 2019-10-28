//
// Trajectory.cpp
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


#include "Trajectory.hpp"
#include "Simulator.hpp"
#include <iostream>


using namespace std;


//
// Trajectory
//


std::string Trajectory::class_name() const
{
    cerr << "[Trajectory] Warning: virtual class_name() has not been defined in derived class.\n";
    return "Trajectory";
}


Parameters Trajectory::parameters() const
{
    cerr << "[Trajectory] Warning: virtual parameters() has not been defined in derived class.\n";
    return Parameters();
}


void Trajectory::configure(const Parameters& parameters, const Registry& registry)
{
    cerr << "[Trajectory] Warning: virtual configure() has not been defined in derived class.\n";
}


// 
// Trajectory_PopulationComposite
//


Trajectory_PopulationComposite::Trajectory_PopulationComposite(const string& id, 
                                                               const TrajectoryPtrs& trajectories)
:   Trajectory(id), trajectories_(trajectories)
{}


double Trajectory_PopulationComposite::value(size_t generation_index,
                                             size_t population_index) const
{
    if (population_index >= trajectories_.size())
        throw runtime_error("[Trajectory_PopulationComposite] Bad population_index.");
    return trajectories_[population_index]->value(generation_index, population_index);
}


Parameters Trajectory_PopulationComposite::parameters() const
{
    Parameters parameters;
    ostringstream trajectories_string;
    for (TrajectoryPtrs::const_iterator it=trajectories_.begin(); it!=trajectories_.end(); ++it)
        trajectories_string << (*it)->object_id() << " ";
    parameters.insert_name_value("trajectories", trajectories_string.str());
    return parameters;    
}


void Trajectory_PopulationComposite::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> trajectory_ids = parameters.value_vector<string>("trajectories");    
    for (vector<string>::const_iterator trajectory_id=trajectory_ids.begin(); trajectory_id!=trajectory_ids.end(); ++trajectory_id)
        trajectories_.push_back(registry.get<Trajectory>(*trajectory_id));
}


void Trajectory_PopulationComposite::initialize(const SimulatorConfig& config)
{
    if (trajectories_.size() != config.population_config_generator->population_count())
        throw runtime_error("[Trajectory_PopulationComposite] Trajectory count does not match population count.");
}


void Trajectory_PopulationComposite::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (TrajectoryPtrs::const_iterator it=trajectories_.begin(); it!=trajectories_.end(); ++it)
    {
        (*it)->write_configuration(os, ids_written);
        ids_written.insert((*it)->object_id());
    }
}


// 
// Trajectory_GenerationComposite
//


Trajectory_GenerationComposite::Trajectory_GenerationComposite(const string& id, 
                                                               const GenerationTrajectoryMap& trajectories)
:   Trajectory(id), trajectories_(trajectories)
{}


double Trajectory_GenerationComposite::value(size_t generation_index,
                                             size_t population_index) const
{
    if (trajectories_.empty())
        throw runtime_error("[Trajectory_GenerationComposite] Initialization error: no trajectories.\n"
                            "Check parameter name: generation:trajectory");

    if (!trajectories_.count(0))
        throw runtime_error("[Trajectory_GenerationComposite] Initialization error: generation 0 not specified.");

    GenerationTrajectoryMap::const_iterator it = trajectories_.upper_bound(generation_index);
    if (it != trajectories_.begin()) --it;

    const Trajectory& trajectory = *it->second;
    
    return trajectory.value(generation_index, population_index);
}


Parameters Trajectory_GenerationComposite::parameters() const
{
    Parameters parameters;
    for (GenerationTrajectoryMap::const_iterator it=trajectories_.begin(); it!=trajectories_.end(); ++it)
        parameters.insert_name_value("generation:trajectory", boost::lexical_cast<string>(it->first) + " " + it->second->object_id());
    return parameters;    
}


void Trajectory_GenerationComposite::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> parameter_values = parameters.values<string>("generation:trajectory");    

    for (vector<string>::const_iterator it=parameter_values.begin(); it!=parameter_values.end(); ++it)
    {
        istringstream iss(*it);
        size_t generation_index;
        string trajectory_id;
        iss >> generation_index >> trajectory_id;
        trajectories_[generation_index] = registry.get<Trajectory>(trajectory_id);
    }
}


void Trajectory_GenerationComposite::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (GenerationTrajectoryMap::const_iterator it=trajectories_.begin(); it!=trajectories_.end(); ++it)
    {
        it->second->write_configuration(os, ids_written);
        ids_written.insert(it->second->object_id());
    }
}


// 
// Trajectory_Constant
//


Parameters Trajectory_Constant::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("value", value_);
    return parameters;
}


void Trajectory_Constant::configure(const Parameters& parameters, const Registry& registry)
{
    value_ = parameters.value<double>("value");
}


// 
// Trajectory_Linear
//


Parameters Trajectory_Linear::parameters() const
{
    Parameters parameters;

    if (generation_index_begin_ != -1ul)
    {
        ostringstream oss_begin;
        oss_begin << generation_index_begin_ << " " << value_begin_;
        ostringstream oss_end;
        oss_end << generation_index_end_ << " " << value_end_;
        parameters.insert_name_value("begin:value", oss_begin.str());
        parameters.insert_name_value("end:value", oss_end.str());
    }
    else
    {
        parameters.insert_name_value("slope", slope_);
        parameters.insert_name_value("intercept", intercept_);
    }

    return parameters;
}


void Trajectory_Linear::configure(const Parameters& parameters, const Registry& registry)
{
    if (parameters.count("begin:value"))
    {
        istringstream iss_begin(parameters.value<string>("begin:value"));
        iss_begin >> generation_index_begin_ >> value_begin_;
        istringstream iss_end(parameters.value<string>("end:value"));
        iss_end >> generation_index_end_ >> value_end_;
        slope_ = (value_end_ - value_begin_) / (generation_index_end_ - generation_index_begin_);
        intercept_ = value_begin_ - slope_ * generation_index_begin_;
    }
    else
    {
        slope_ = parameters.value<double>("slope");
        intercept_ = parameters.value<double>("intercept");
    }
}


//
// Trajectory_Exponential
//


double Trajectory_Exponential::value(size_t generation_index, size_t population_index) const
{
    if (generation_index < generation_index_begin_)
    {
        ostringstream oss;
        oss << "[Trajectory_Exponential] Undefined at generation " << generation_index;
        throw runtime_error(oss.str().c_str());
    }

    return value_begin_ * exp(rate_ * (generation_index - generation_index_begin_));
}


Parameters Trajectory_Exponential::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("generation_begin", generation_index_begin_);
    parameters.insert_name_value("value_begin", value_begin_);
    parameters.insert_name_value("rate", rate_);
    return parameters;
}


void Trajectory_Exponential::configure(const Parameters& parameters, const Registry& registry)
{
    generation_index_begin_ = parameters.value<size_t>("generation_begin");
    value_begin_ = parameters.value<double>("value_begin");
    rate_ = parameters.value<double>("rate");
}


