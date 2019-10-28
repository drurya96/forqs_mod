//
// RecombinationPositionGeneratorImplementation.cpp
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


#include "RecombinationPositionGeneratorImplementation.hpp"
#include "Simulator.hpp"
#include <stdexcept>


using namespace std;


//
// RecombinationPositionGenerator_Trivial
//


vector<unsigned int> RecombinationPositionGenerator_Trivial::get_positions(size_t chromosome_pair_index) const
{
    vector<unsigned int> result;
    if (Random::uniform_01()>=.5) result.push_back(0); // start with 2nd chromosome
    return result;
}


Parameters RecombinationPositionGenerator_Trivial::parameters() const
{
    return Parameters();
}


void RecombinationPositionGenerator_Trivial::configure(const Parameters& parameters, const Registry& registry)
{}


//
// RecombinationPositionGenerator_SingleCrossover
//


vector<unsigned int> RecombinationPositionGenerator_SingleCrossover::get_positions(size_t chromosome_pair_index) const
{
    if (chromosome_pair_index >= chromosome_lengths_.size())
        throw runtime_error("[RecombinationPositionGenerator_SingleCrossover] Invalid chromosome_pair_index.");

    unsigned int chromosome_length = chromosome_lengths_[chromosome_pair_index];

    vector<unsigned int> result;
    if (Random::uniform_01()>=.5) result.push_back(0); // start with 2nd chromosome

    if (Random::uniform_01()>=.5) return result; // no recombination

    result.push_back(Random::uniform_integer(0, chromosome_length));

    return result;
}


Parameters RecombinationPositionGenerator_SingleCrossover::parameters() const
{
    return Parameters();
}


void RecombinationPositionGenerator_SingleCrossover::configure(const Parameters& parameters, const Registry& registry)
{}


void RecombinationPositionGenerator_SingleCrossover::initialize(const SimulatorConfig& config)
{
    chromosome_lengths_ = config.population_config_generator->chromosome_lengths();
}


//
// RecombinationPositionGenerator_Uniform
//


RecombinationPositionGenerator_Uniform::ChromosomeInfo::ChromosomeInfo(size_t _length, double _rate)
:   length(_length), rate(_rate)
{
    poisson = Random::create_poisson_distribution("", rate);
}


RecombinationPositionGenerator_Uniform::RecombinationPositionGenerator_Uniform(const string& id, 
                                                                               vector<ChromosomeInfo> infos)
:   RecombinationPositionGenerator(id), 
    common_rate_(1.0),
    infos_(infos)
{}


vector<unsigned int> RecombinationPositionGenerator_Uniform::get_positions(size_t chromosome_pair_index) const
{
    if (chromosome_pair_index >= infos_.size())
        throw runtime_error("[RecombinationPositionGenerator_Uniform] Invalid chromosome_pair_index.");

    const ChromosomeInfo& info = infos_[chromosome_pair_index];

    size_t recombination_count = (size_t)info.poisson->random_value();
    vector<unsigned int> result;
    result.reserve(recombination_count+1);

    if (Random::uniform_01()>=.5) result.push_back(0); // start with 2nd chromosome

    vector<size_t> positions = Random::random_indices_without_replacement(info.length, recombination_count);

    for (vector<size_t>::const_iterator position=positions.begin(); position!=positions.end(); ++position)
        if (*position) result.push_back(*position);

    return result;
}


Parameters RecombinationPositionGenerator_Uniform::parameters() const
{
    bool rates_equal = true;

    for (ChromosomeInfos::const_iterator info=infos_.begin(); info!=infos_.end(); ++info)
    {
        if (info->rate != common_rate_)
        {
            rates_equal = false;
            break;
        }
    }

    Parameters parameters;

    if (rates_equal)
        parameters.insert_name_value("rate", common_rate_);
    else
    {
        ostringstream rates;
        for (ChromosomeInfos::const_iterator info=infos_.begin(); info!=infos_.end(); ++info)
            rates << info->rate << " ";
        parameters.insert_name_value("rates", rates.str());
    }

    return parameters;
}


void RecombinationPositionGenerator_Uniform::configure(const Parameters& parameters, const Registry& registry)
{
    common_rate_ = parameters.value<double>("rate", 1.0);

    if (parameters.count("rates"))
    {
        vector<double> rates = parameters.value_vector<double>("rates");

        for (vector<double>::const_iterator it=rates.begin(); it!=rates.end(); ++it)
            infos_.push_back(ChromosomeInfo(0, *it)); // we don't know the chromosome lengths until initialize()
    }
}


void RecombinationPositionGenerator_Uniform::initialize(const SimulatorConfig& config)
{
    // configure() is called before initialize(); we have rates, but no chromosome lengths

    const size_t chromosome_pair_count = config.population_config_generator->chromosome_pair_count();
    const vector<unsigned int>& chromosome_lengths = config.population_config_generator->chromosome_lengths();

    if (chromosome_pair_count != chromosome_lengths.size())
    {
        ostringstream oss;
        oss << "[RecombinationPositionGenerator_Uniform] Number of chromosome_lengths " 
            << "doesn't match chromosome_pair_count in "
            << config.population_config_generator->class_name() << " "
            << config.population_config_generator->object_id();
        throw runtime_error(oss.str().c_str());
    }

    if (infos_.empty())
        for (size_t i=0; i<chromosome_lengths.size(); ++i)
            infos_.push_back(ChromosomeInfo(0, common_rate_));

    if (infos_.size() != chromosome_lengths.size())
        throw runtime_error("[RecombinationPositionGenerator_Uniform] Number of rates doesn't match number of chromosomes.");

    vector<unsigned int>::const_iterator chromosome_length = chromosome_lengths.begin();
    for (ChromosomeInfos::iterator info=infos_.begin(); info!=infos_.end(); ++info, ++chromosome_length)
    {
        if (*chromosome_length == 0)
            throw runtime_error("[RecombinationPositionGenerator_Uniform] Chromosome length must be nonzero.");

        info->length = *chromosome_length;
    }

	//cout << "Initiliazing a thing..." << endl;

}


//
// RecombinationPositionGenerator_RecombinationMap
//


vector<unsigned int> RecombinationPositionGenerator_RecombinationMap::get_positions(size_t chromosome_pair_index) const
{
    if (chromosome_pair_index >= recombination_maps_.size())
        throw runtime_error("[RecombinationPositionGenerator_RecombinationMap::get_positions()] Index out of bounds.");

    vector<unsigned int> positions = recombination_maps_[chromosome_pair_index]->random_positions(); // may not be sorted 
    if (Random::uniform_01()>=.5) positions.push_back(0); // start with 2nd chromosome
    sort(positions.begin(), positions.end());

    return positions; 
}


RecombinationPositionGenerator_RecombinationMap::RecombinationPositionGenerator_RecombinationMap(
    const string& id,
    const vector<string>& filenames)
:   RecombinationPositionGenerator(id), filenames_(filenames)
{
    if (!filenames_.empty())
        read_files();
}


void RecombinationPositionGenerator_RecombinationMap::read_files()
{
    recombination_maps_.clear();

    for (vector<string>::const_iterator it=filenames_.begin(); it!=filenames_.end(); ++it)
        recombination_maps_.push_back(boost::shared_ptr<RecombinationMap>(
            new RecombinationMap(*it)));
}


Parameters RecombinationPositionGenerator_RecombinationMap::parameters() const
{
    Parameters parameters;
    for (vector<string>::const_iterator it=filenames_.begin(); it!=filenames_.end(); ++it)
        parameters.insert_name_value("filename", *it);
    return parameters;
}


void RecombinationPositionGenerator_RecombinationMap::configure(const Parameters& parameters, const Registry& registry)
{
    filenames_ = parameters.values<string>("filename");
    read_files();
}


//
// RecombinationPositionGenerator_Composite
//


RecombinationPositionGenerator_Composite::RecombinationPositionGenerator_Composite(
    const string& id)
:   RecombinationPositionGenerator(id),
    default_rpg_(new RecombinationPositionGenerator_Trivial("id_rpg_composite_default_trivial"))
{}


vector<unsigned int> RecombinationPositionGenerator_Composite::get_positions(size_t chromosome_pair_index) const
{
    if (rpg_map_.count(chromosome_pair_index))
        return rpg_map_.at(chromosome_pair_index)->get_positions(chromosome_pair_index);

    return default_rpg_->get_positions(chromosome_pair_index);
}


Parameters RecombinationPositionGenerator_Composite::parameters() const
{
    Parameters parameters;

    parameters.insert_name_value("default_recombination_position_generator", 
        default_rpg_->object_id());

    for (RPGMap::const_iterator it=rpg_map_.begin(); it!=rpg_map_.end(); ++it)
    {
        ostringstream value;
        value << it->first + 1 << " " << it->second->object_id(); // 1-based chromosome parameter
        parameters.insert_name_value("chromosome:recombination_position_generator", value.str());
    }

    return parameters;
}


void RecombinationPositionGenerator_Composite::configure(const Parameters& parameters, const Registry& registry)
{
    if (parameters.count("default_recombination_position_generator"))
        default_rpg_ = registry.get<RecombinationPositionGenerator>(
            parameters.value<string>("default_recombination_position_generator"));

    vector<string> configurations = parameters.values<string>("chromosome:recombination_position_generator");
    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        istringstream iss(*it);
        size_t chromosome = 0;
        string id;
        iss >> chromosome >> id;

        if (chromosome == 0)
            throw runtime_error("[RecombinationPositionGenerator_Composite] Invalid chromosome 0.");
        size_t chromosome_pair_index = chromosome - 1; // 0-based index

        if (rpg_map_.count(chromosome_pair_index))
            throw runtime_error("[RecombinationPositionGenerator_Composite] Duplicate chromosome specification.");

        RecombinationPositionGeneratorPtr rpg = registry.get<RecombinationPositionGenerator>(id);

        rpg_map_[chromosome_pair_index] = rpg;
    }
}


