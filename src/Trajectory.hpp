//
// Trajectory.hpp
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


#ifndef _TRAJECTORY_HPP_
#define _TRAJECTORY_HPP_


#include "Configurable.hpp"
#include "shared_ptr.hpp"


///
/// \defgroup Trajectories Trajectories
///
/// classes for defining values that are variable across populations or generations,
/// i.e. real-valued functions of (generation_index, population_index)
///


//
// Trajectory
//

///
/// interface class
///

class Trajectory : public Configurable
{
    public:

    virtual double value(size_t generation_index, size_t population_index) const = 0;

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    protected:

    Trajectory(const std::string& id) : Configurable(id) {}
};


typedef shared_ptr<Trajectory> TrajectoryPtr;
typedef std::vector<TrajectoryPtr> TrajectoryPtrs;


//
// Trajectory_PopulationComposite
//

///
/// per-population specification of trajectories
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// trajectories = \<id_trajectory_1\> [\<id_trajectory_2\> ...] | none | trajectory count must match population count
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup Trajectories
///


class Trajectory_PopulationComposite : public Trajectory
{
    public:

    Trajectory_PopulationComposite(const std::string& id, 
                                   const TrajectoryPtrs& trajectories = TrajectoryPtrs()); 

    virtual double value(size_t generation_index, size_t population_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "Trajectory_PopulationComposite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    TrajectoryPtrs trajectories_;
};


//
// Trajectory_GenerationComposite
//

///
/// per-generation specification of trajectories
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// generation:trajectory = \<int_generation_begin\> \<id_trajectory\> | none | multiple allowed; generation 0 required
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup Trajectories
///


class Trajectory_GenerationComposite : public Trajectory
{
    public:

    typedef std::map<size_t, TrajectoryPtr> GenerationTrajectoryMap;

    Trajectory_GenerationComposite(const std::string& id, 
                                   const GenerationTrajectoryMap& trajectories = GenerationTrajectoryMap()); 

    virtual double value(size_t generation_index, size_t population_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "Trajectory_GenerationComposite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    GenerationTrajectoryMap trajectories_; // generation_index -> Trajectory
};


//
// Trajectory_Constant
//

///
/// trajectory with a constant value
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// value = \<float\> | none | required
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup Trajectories
///


class Trajectory_Constant : public Trajectory
{
    public:

    Trajectory_Constant(const std::string& id, double value = 0.)
    :   Trajectory(id), value_(value)
    {}

    virtual double value(size_t generation_index, size_t population_index) const {return value_;}

    // Configurable interface

    virtual std::string class_name() const {return "Trajectory_Constant";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:
    double value_;
};


//
// Trajectory_Linear
//

///
/// trajectory with values linear in generation_index 
/// 
/// 2 parametrizations are available:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// slope = \<float\> | none | required
/// intercept = \<float\> | none | required
///
/// or:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// begin:value = \<int_generation_begin\> <float_value\> | none | required
/// end:value = \<int_generation_end\> <float_value\> | none | required
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup Trajectories
///


class Trajectory_Linear : public Trajectory
{
    public:

    Trajectory_Linear(const std::string& id, double slope = 0., double intercept = 0.)
    :   Trajectory(id), slope_(slope), intercept_(intercept),
        generation_index_begin_(-1ul), generation_index_end_(-1ul), 
        value_begin_(0), value_end_(0)
    {}

    virtual double value(size_t generation_index, size_t population_index) const 
    {
        return slope_ * generation_index + intercept_;
    }

    // Configurable interface

    virtual std::string class_name() const {return "Trajectory_Linear";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    double slope_;
    double intercept_;

    size_t generation_index_begin_;
    size_t generation_index_end_;
    double value_begin_;
    double value_end_;
};


//
// Trajectory_Exponential
//

///
/// trajectory with values exponential in generation
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// generation_begin = \<float\> | none | required
/// value_begin = \<float\> | none | required
/// rate = \<float\> | none | required
///
/// If t is the generation, returns value_begin * exp[rate * (t - generation_begin)]
///
/// Example: [example_trajectories.txt](../../examples/example_trajectories.txt)
///
/// \ingroup Trajectories
///


class Trajectory_Exponential : public Trajectory
{
    public:

    Trajectory_Exponential(const std::string& id, 
                           size_t generation_index_begin = 0,
                           double value_begin = 0., 
                           double rate = 0.)
    :   Trajectory(id),
        generation_index_begin_(generation_index_begin),
        value_begin_(value_begin),
        rate_(rate)
    {}

    virtual double value(size_t generation_index, size_t population_index) const;

    // Configurable interface

    virtual std::string class_name() const {return "Trajectory_Exponential";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    size_t generation_index_begin_;
    double value_begin_;
    double rate_;
};


#endif // _TRAJECTORY_HPP_

