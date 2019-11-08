//
// FitnessFunctionImplementation.hpp
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


#ifndef _FITNESSFUNCTIONIMPLEMENTATION_HPP_
#define _FITNESSFUNCTIONIMPLEMENTATION_HPP_


#include "QuantitativeTrait.hpp"
#include "PopulationData.hpp"
#include "shared_ptr.hpp"
#include <stdexcept>


///
/// \defgroup FitnessFunctions FitnessFunctions
///
/// classes for calculating fitnesses from quantitative traits
///


//
// FitnessFunction_Trivial
//

///
/// returns fitness 1 for all individuals
/// 
/// Parameters: none
///
/// Example: [example_stepping_stone.txt](../../examples/example_stepping_stone.txt)
///
/// \ingroup FitnessFunctions
///

class FitnessFunction_Trivial : public QuantitativeTrait
{
    public:

    FitnessFunction_Trivial(const std::string& id) : QuantitativeTrait(id) {}

    virtual void calculate_trait_values(const PopulationData& population_data) const
    {
        (*population_data.trait_values)[object_id()] = 
            DataVectorPtr(new DataVector(population_data.population_size, 1));
    }

    // Configurable interface

    virtual std::string class_name() const {return "FitnessFunction_Trivial";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}
};


//
// FitnessFunction_Optimum
//

///
/// Fitness function assigning fitness 1 to optimum trait value, with fitness 
/// decaying away from the optimum.
///
/// The decay can be polynomial:
/// - fitness == (1 - abs(trait_value - optimum)/radius)^power
///
/// or Gaussian:
/// - fitness = exp[ - (trait_value - optimum)^2 / 2*width^2 ]
///
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_trait = \<id\> | none | required
/// optimum = \<double\> | none | required
/// radius:power = \<double\> \<double\> | none | optional*
/// gaussian_width = \<double\> | none | optional*
///
/// *Exactly one of radius:power or gaussian_width must be specified.
///
/// Example: [example_stepping_stone.txt](../../examples/example_stepping_stone.txt)
///
/// Note: the term 'width' to parametrize Gaussian fitness decay comes from
///       Lande 1976 Evolution (Natural Selection and Random Genetic Drift in Phenotypic Evolution)
///
/// \ingroup FitnessFunctions
///

class FitnessFunction_Optimum : public QuantitativeTrait
{
    public:

    FitnessFunction_Optimum(const std::string& id, const std::string& quantitative_trait_id = "",
                            double optimum = 0, double radius = 0, double power = 0);

    virtual void calculate_trait_values(const PopulationData& population_data) const;

    // Configurable interface

    virtual std::string class_name() const {return "FitnessFunction_Optimum";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::string qtid_;
    double optimum_;
    double radius_;
    double power_;
    double gaussian_width_;
};


//
// FitnessFunction_TruncationSelection
//

///
/// Fitness function assigning fitness 1 if trait value exceeds threshold, 0 otherwise.
/// Threshold is determined each generation by the user-specified proportion of 
/// individuals to be selected to reproduce.
///
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_trait = \<id\> | none | required
/// proportion_selected = \<float\> | none | required
/// lower_tail = 1 | 0 | select from the lower tail of the distribution
/// single_threshold_population = \<int\> | none | specify population for single threshold
/// ignore_zero = 1 | 0 | select proportion of non-zero values only
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup FitnessFunctions
///

class FitnessFunction_TruncationSelection : public QuantitativeTrait
{
    public:

    FitnessFunction_TruncationSelection(const std::string& id);

    void calculate_trait_values_with_threshold(const PopulationData& population_data, double threshold) const;

    virtual void calculate_trait_values(const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "FitnessFunction_TruncationSelection";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::string qtid_;
    double proportion_selected_;
    bool lower_tail_;
    bool single_threshold_;
    size_t single_threshold_population_index_;
    bool ignore_zero_;

    double calculate_threshold(const PopulationData& population_data) const;
};

// Austin Drury 10/8/2019
//
// FitnessFunction_BoundedSelection
//

///
/// Fitness function assigning fitness 1 if trait value falls between boundaries, 0 otherwise.
/// Boundaries are set by the user, either upper bound or lower bound.
///
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_trait = \<id\> | none | required
/// proportion_selected = \<float\> | none | required
/// lower_tail = 1 | 0 | select from the lower tail of the distribution
/// single_threshold_population = \<int\> | none | specify population for single threshold
/// ignore_zero = 1 | 0 | select proportion of non-zero values only
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup FitnessFunctions
///

class FitnessFunction_BoundedSelection : public QuantitativeTrait
{
    public:

    FitnessFunction_BoundedSelection(const std::string& id);

    void calculate_trait_values_with_threshold(const PopulationData& population_data) const;

    virtual void calculate_trait_values(const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "FitnessFunction_BoundedSelection";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::string qtid_;
	double lower_bound_;
	double upper_bound_;
};

// Austin Drury 10/15/2019
//
// FitnessFunction_Recombination
//

///
/// Fitness function assigning fitness based on distributions around maximum and minimum values
/// Boundaries are set by the user, either upper bound or lower bound or both. User also sets varation
///
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_trait = \<id\> | none | required
/// proportion_selected = \<float\> | none | required
/// lower_tail = 1 | 0 | select from the lower tail of the distribution
/// single_threshold_population = \<int\> | none | specify population for single threshold
/// ignore_zero = 1 | 0 | select proportion of non-zero values only
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup FitnessFunctions
///

class FitnessFunction_Recombination : public QuantitativeTrait
{
    public:

    FitnessFunction_Recombination(const std::string& id);

	//void set_recombination_array();

    void modify_trait_values(const PopulationData& population_data) const;

	double calculate_new_value(const double value) const;

    virtual void calculate_trait_values(const PopulationDataPtrs& population_datas) const;

    // Configurable interface

    virtual std::string class_name() const {return "FitnessFunction_Recombination";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

	double getLowerBound(void){ return this->lower_bound_; }
	double getUpperBound(void){ return this->upper_bound_; }

    //private:

    std::string qtid_;
	double lower_bound_;
	double upper_bound_;
	double variation_;
};



#endif // _FITNESSFUNCTIONIMPLEMENTATION_HPP_

