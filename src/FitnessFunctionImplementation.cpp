//
// FitnessFunctionImplementation.cpp
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


#include "FitnessFunctionImplementation.hpp"
#include "Simulator.hpp"


using namespace std;


//
// FitnessFunction_Optimum
//


namespace {
const double invalid_gaussian_width_ = -1.0;
} // namespace


FitnessFunction_Optimum::FitnessFunction_Optimum(const string& id, 
                                                 const string& quantitative_trait_id, 
                                                 double optimum, 
                                                 double radius, 
                                                 double power)
:   QuantitativeTrait(id), qtid_(quantitative_trait_id),
    optimum_(optimum), radius_(radius), power_(power),
    gaussian_width_(invalid_gaussian_width_)
{}


void FitnessFunction_Optimum::calculate_trait_values(const PopulationData& population_data) const
{
    const DataVector& trait_values = *population_data.trait_values->get(qtid_);
    DataVectorPtr fitnesses(new DataVector(trait_values.size()));

    DataVector::const_iterator trait_value = trait_values.begin();
    DataVector::iterator fitness = fitnesses->begin();

    if (gaussian_width_ == invalid_gaussian_width_)
    {
        for (; trait_value!=trait_values.end(); ++trait_value, ++fitness)
        {
            // polynomial fitness == (1 - |trait_value - optimum|/radius)^power

            double d = min(fabs(*trait_value - optimum_)/radius_, 1.0); // d == scaled distance from optimum in [0,1]
            *fitness = pow((1 - d), power_);
        }
    }
    else
    {
        // Gaussian fitness = exp[ - (trait_value - optimum)^2 / 2*width^2 ]

        for (; trait_value!=trait_values.end(); ++trait_value, ++fitness)
        {
            double d = (*trait_value - optimum_)/gaussian_width_;
            *fitness = exp(-.5 * pow(d,2.));
        }
    }

    (*population_data.trait_values)[object_id()] = fitnesses;
}


Parameters FitnessFunction_Optimum::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid_);
    parameters.insert_name_value("optimum", optimum_);

    if (gaussian_width_ == invalid_gaussian_width_)
    {
        ostringstream radius_power;
        radius_power << radius_ << " " << power_;
        parameters.insert_name_value("radius:power", radius_power.str());
    }
    else
    {
        parameters.insert_name_value("gaussian_width", gaussian_width_);
    }

    return parameters;
}


void FitnessFunction_Optimum::configure(const Parameters& parameters, const Registry& registry)
{
    if (!parameters.count("gaussian_width") && !parameters.count("radius:power") ||
        parameters.count("gaussian_width") && parameters.count("radius:power"))
        throw runtime_error("[FitnessFunction_Optimum] Exactly one of radius:power or gaussian_width must be specified.");

    qtid_ = parameters.value<string>("quantitative_trait");
    optimum_ = parameters.value<double>("optimum");

    if (parameters.count("gaussian_width"))
    {
        gaussian_width_ = parameters.value<double>("gaussian_width");
    }
    else
    {
        vector<double> radius_power = parameters.value_vector<double>("radius:power");
        if (radius_power.size() != 2)
            throw runtime_error("[FitnessFunction_Optimum] Invalid specification of radius:power.");

        radius_ = radius_power[0];
        power_ = radius_power[1];
    }
}


// ======================================================================================================================================
// FitnessFunction_TruncationSelection
// ======================================================================================================================================


FitnessFunction_TruncationSelection::FitnessFunction_TruncationSelection(const string& id)
:   QuantitativeTrait(id), 
    proportion_selected_(0),
    lower_tail_(false),
    single_threshold_(false),
    single_threshold_population_index_(0),
    ignore_zero_(false)
{}


void FitnessFunction_TruncationSelection::calculate_trait_values_with_threshold(const PopulationData& population_data, 
                                                                                double threshold) const
{
    const DataVector& trait_values = *population_data.trait_values->get(qtid_);

    // fitness = (trait_value >= threshold) ? 1 : 0;

    DataVectorPtr fitnesses(new DataVector(trait_values.size()));

    if (lower_tail_)
    {
        if (!ignore_zero_)
        {
            transform(trait_values.begin(), trait_values.end(), fitnesses->begin(),
                      bind2nd(less_equal<double>(),threshold));
        }
        else
        {
            DataVector::const_iterator it=trait_values.begin();
            for (DataVector::iterator jt=fitnesses->begin(); jt!=fitnesses->end(); ++jt, ++it)
                *jt = (*it <= threshold && *it != 0);
        }
    }
    else
    {
        if (!ignore_zero_)
        {
            transform(trait_values.begin(), trait_values.end(), fitnesses->begin(),
                      bind2nd(greater_equal<double>(),threshold));
        }
        else
        {
            // note: there is no unit test for this case: negative trait values, upper tail, ignore_zero=1
            DataVector::const_iterator it=trait_values.begin();
            for (DataVector::iterator jt=fitnesses->begin(); jt!=fitnesses->end(); ++jt, ++it)
                *jt = (*it >= threshold && *it != 0);
        }
    }

    // cout << "threshold: " << threshold << endl
    //      << "trait_values: " << trait_values << endl
    //      << "fitnesses: " << *fitnesses << endl;

    (*population_data.trait_values)[object_id()] = fitnesses;
}


void FitnessFunction_TruncationSelection::calculate_trait_values(const PopulationDataPtrs& population_datas) const
{
    if (single_threshold_ && single_threshold_population_index_ >= population_datas.size())
    {
        QuantitativeTraitPtr qt(new FitnessFunction_Trivial(object_id()));
        qt->calculate_trait_values(population_datas);
        return;
    }

    vector<double> thresholds;

    if (single_threshold_)
    {
        double threshold = calculate_threshold(*population_datas[single_threshold_population_index_]);
        thresholds = vector<double>(population_datas.size(), threshold);
    }
    else // separate threshold for each population
    {
        for (PopulationDataPtrs::const_iterator popdata = population_datas.begin();
             popdata != population_datas.end(); ++popdata)
            thresholds.push_back(calculate_threshold(**popdata));
    }

    vector<double>::const_iterator threshold = thresholds.begin(); 
    for (PopulationDataPtrs::const_iterator popdata = population_datas.begin(); 
         popdata != population_datas.end(); ++popdata, ++threshold)
    {
        if ((*popdata)->trait_values->count(object_id()))
            throw runtime_error("[QuantitativeTrait::calculate_trait_values()] Quantitative trait id already used.");
        calculate_trait_values_with_threshold(**popdata, *threshold);
    }
}


Parameters FitnessFunction_TruncationSelection::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid_);
    parameters.insert_name_value("proportion_selected", proportion_selected_);
    if (lower_tail_) parameters.insert_name_value("lower_tail", lower_tail_);
    if (single_threshold_) 
        parameters.insert_name_value("single_threshold_population", single_threshold_population_index_ + 1); // 1-based
    if (ignore_zero_) parameters.insert_name_value("ignore_zero", ignore_zero_);
    return parameters;
}

void FitnessFunction_TruncationSelection::configure(const Parameters& parameters, const Registry& registry)
{
    qtid_ = parameters.value<string>("quantitative_trait");
    proportion_selected_ = parameters.value<double>("proportion_selected");

    lower_tail_ = parameters.value<bool>("lower_tail", false);

    single_threshold_ = parameters.count("single_threshold_population");
    if (single_threshold_)
    {
        int population = parameters.value<int>("single_threshold_population"); // 1-based
        if (population < 1)
            throw runtime_error("[FitnessFunction_TruncationSelection] Invalid population." );
        single_threshold_population_index_ = size_t(population - 1); // 0-based
    }

    ignore_zero_ = parameters.value<bool>("ignore_zero", false);
}


double FitnessFunction_TruncationSelection::calculate_threshold(const PopulationData& population_data) const
{
    const DataVector& trait_values = *population_data.trait_values->get(qtid_);

    // copy the trait values to sort and obtain threshold

    DataVector trait_values_copy;

    if (ignore_zero_)
    {
        remove_copy_if(trait_values.begin(), trait_values.end(), back_inserter(trait_values_copy),
                       bind2nd(equal_to<double>(), 0.));
    }
    else
    {
        trait_values_copy.reserve(trait_values.size());
        copy(trait_values.begin(), trait_values.end(), back_inserter(trait_values_copy));
    }        

    if (trait_values_copy.empty())
        throw runtime_error("[FitnessFunction_TruncationSelection " + object_id() + "] Empty trait values.");

    size_t index_cutoff = static_cast<size_t>(floor(trait_values_copy.size() * proportion_selected_));
    if (index_cutoff > 0) index_cutoff -= 1;

    // linear time partial sort puts correct element in position index_cutoff

    if (lower_tail_)
        nth_element(trait_values_copy.begin(), trait_values_copy.begin() + index_cutoff,
                    trait_values_copy.end(), less<double>());
    else
        nth_element(trait_values_copy.begin(), trait_values_copy.begin() + index_cutoff,
                    trait_values_copy.end(), greater<double>());

    const double threshold = trait_values_copy[index_cutoff];

    return threshold;
}

// ======================================================================================================================================
// FitnessFunction_BoundedSelection
// Created by Austin Drury 10/19
// ======================================================================================================================================

// default lower bound is zero (no negative fitness)
// default upper bound is 1000000000


// calculate_threshold above returns the absolute fitness to accept.
// This function will be directly given an absolute fitness so calculate threshold is not needed

FitnessFunction_BoundedSelection::FitnessFunction_BoundedSelection(const string& id)
:   QuantitativeTrait(id), 
    lower_bound_(0),
	upper_bound_(1000000000)
    //lower_tail_(false),
    //single_threshold_(false),
    //single_threshold_population_index_(0),
    //ignore_zero_(false)
{}


void FitnessFunction_BoundedSelection::calculate_trait_values_with_threshold(const PopulationData& population_data) const
{
    const DataVector& trait_values = *population_data.trait_values->get(qtid_);

    // fitness = (trait_value >= threshold) ? 1 : 0;

    DataVectorPtr fitnesses_lower(new DataVector(trait_values.size()));
	DataVectorPtr fitnesses_upper(new DataVector(trait_values.size()));
	DataVectorPtr fitnesses(new DataVector(trait_values.size()));

    transform(trait_values.begin(), trait_values.end(), fitnesses_upper->begin(), 
		bind2nd(less_equal<double>(),upper_bound_));
	transform(trait_values.begin(), trait_values.end(), fitnesses_lower->begin(), 
		bind2nd(greater_equal<double>(),lower_bound_));

	vector<double>::const_iterator lt = fitnesses_lower->begin();
	vector<double>::iterator ft = fitnesses->begin();

	for (vector<double>::const_iterator ut = fitnesses_upper->begin(); ut != fitnesses_upper->end(); ++ut, ++lt, ++ft){
		if(*lt && *ut){
			*ft = 1;
		}
		else{
			*ft = 0;
		}
	} 

    // cout << "threshold: " << threshold << endl
    //      << "trait_values: " << trait_values << endl
    //      << "fitnesses: " << *fitnesses << endl;

    (*population_data.trait_values)[object_id()] = fitnesses;
}


void FitnessFunction_BoundedSelection::calculate_trait_values(const PopulationDataPtrs& population_datas) const
{

    for (PopulationDataPtrs::const_iterator popdata = population_datas.begin(); 
         popdata != population_datas.end(); ++popdata)
    {
        if ((*popdata)->trait_values->count(object_id()))
            throw runtime_error("[QuantitativeTrait::calculate_trait_values()] Quantitative trait id already used.");
        calculate_trait_values_with_threshold(**popdata);
    }
}


Parameters FitnessFunction_BoundedSelection::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid_);
    parameters.insert_name_value("lower_bound", lower_bound_);
	parameters.insert_name_value("upper_bound", upper_bound_);
    return parameters;
}


// For now assume lower bound is >= 0 AND upper bount is <= 1000000000 AND upper bound > lower bound
void FitnessFunction_BoundedSelection::configure(const Parameters& parameters, const Registry& registry)
{
    qtid_ = parameters.value<string>("quantitative_trait");
    lower_bound_ = parameters.value<double>("lower_bound", 0);
	upper_bound_ = parameters.value<double>("upper_bound", 1000000000);
}



// ======================================================================================================================================
// FitnessFunction_Recombination
// Created by Austin Drury 10/19
// ======================================================================================================================================

// default lower bound is zero (no negative fitness)
// default upper bound is 1000000000
// default varation is 0.1

FitnessFunction_Recombination::FitnessFunction_Recombination(const string& id)
:   QuantitativeTrait(id), 
    lower_bound_(0),
	upper_bound_(1000000000),
	variation_(0.1)
{}


void FitnessFunction_Recombination::modify_trait_values(const PopulationData& population_data) const
{
    const DataVector& trait_values = *population_data.trait_values->get(qtid_);

    // fitness = (trait_value >= threshold) ? 1 : 0;

	// modify fitness to follow distribution of recombination
	// lower bound function: ((x/lower_bound)^(1/variation))/(1+(x/lower_bound)^(1/varation))
	// upper bound function: 1/(1+(x/lower_bound)^(1/varation))
	//
	// take the minimum fitness?

	double new_value;
	DataVectorPtr fitnesses(new DataVector(trait_values.size()));
	vector<double>::iterator ft = fitnesses->begin();

	for (vector<double>::const_iterator it = trait_values.begin(); it != trait_values.end(); ++it, ++ft){
		//lb = lower_bound_function(*it);
		//ub = upper_bound_function(*it);
		new_value = calculate_new_value(*it);
		*ft = new_value;
	}

    // cout << "threshold: " << threshold << endl
    //      << "trait_values: " << trait_values << endl
    //      << "fitnesses: " << *fitnesses << endl;

    (*population_data.trait_values)[object_id()] = fitnesses;
}

double FitnessFunction_Recombination::calculate_new_value(const double value) const{
	double inside_l = ((double)value-lower_bound_)/variation_;
	double inside_u = ((double)value-upper_bound_)/variation_;
	double lower_v = 0.5 + (tanh(inside_l)*0.5);
	double upper_v = 0.5 + (tanh(inside_u*-1)*0.5);
	//double al = pow((value/lower_bound_),((double)1/variation_));
	//double au = pow((value/upper_bound_),((double)1/variation_));

	//double lower_v = al/((double)1+al);
	//double upper_v = 1/((double)1+au);

	return min(lower_v, upper_v);

}

void FitnessFunction_Recombination::calculate_trait_values(const PopulationDataPtrs& population_datas) const
{

    for (PopulationDataPtrs::const_iterator popdata = population_datas.begin(); 
         popdata != population_datas.end(); ++popdata)
    {
        if ((*popdata)->trait_values->count(object_id()))
            throw runtime_error("[QuantitativeTrait::calculate_trait_values()] Quantitative trait id already used.");
        modify_trait_values(**popdata);
    }
}


Parameters FitnessFunction_Recombination::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("quantitative_trait", qtid_);
    parameters.insert_name_value("lower_bound", lower_bound_);
	parameters.insert_name_value("upper_bound", upper_bound_);
	parameters.insert_name_value("variation", variation_);
    return parameters;
}


// For now assume lower bound is >= 0 AND upper bount is <= 1000000000 AND upper bound > lower bound AND 0 > variation > 1
void FitnessFunction_Recombination::configure(const Parameters& parameters, const Registry& registry)
{
    qtid_ = parameters.value<string>("quantitative_trait");
    lower_bound_ = parameters.value<double>("lower_bound", 0);
	upper_bound_ = parameters.value<double>("upper_bound", 1000000000);
	variation_ = parameters.value<double>("variation", 0.1);
	// This only appears once, so this is a good place to modify the RPGPtrsArray
	cout << "Configuring Fitness Function" << endl;
}

