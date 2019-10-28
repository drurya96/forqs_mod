//
// MutationGeneratorImplementation.cpp
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


#include "MutationGeneratorImplementation.hpp"
#include "Random.hpp"


using namespace std;


//
// MutationGenerator_SingleLocus
//


MutationGenerator_SingleLocus::MutationInfos MutationGenerator_SingleLocus::generate_mutations(
    const Population& population,
    size_t generation_index,
    size_t population_index) const
{
    // note: rate must be calculated each generation, since population size may change
    const double rate = 2 * population.population_size() * mu_; // per chromosome
    Random::DistributionPtr poisson = Random::create_poisson_distribution("id_dummy", rate);

    size_t mutant_count = size_t(poisson->random_value());

    vector<size_t> chromosome_indices = 
        Random::random_indices_without_replacement(2 * population.population_size(), mutant_count);

    MutationInfos result(mutant_count);

    vector<size_t>::const_iterator chromosome_index = chromosome_indices.begin();
    for (MutationInfos::iterator info=result.begin(); info!=result.end(); ++info, ++chromosome_index)
    {
        // individual_index == chromosome_index / 2 
        // which == chromsome_index % 2

        info->individual_index = *chromosome_index >> 1;
        info->locus = locus_;
        info->which = *chromosome_index & 0x01;
        info->value = 1;
    }

    return result;
}


Parameters MutationGenerator_SingleLocus::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("mu", mu_);
    parameters.insert_name_value("locus", locus_.object_id());
    return parameters;
}


void MutationGenerator_SingleLocus::configure(const Parameters& parameters, const Registry& registry)
{
    locus_ = *registry.get<Locus>(parameters.value<string>("locus"));
    mu_ = parameters.value<double>("mu");
}


void MutationGenerator_SingleLocus::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    locus_.write_configuration(os, ids_written);
}


//
// MutationGenerator_Regions
//


string MutationGenerator_Regions::RegionInfo::configuration() const
{
    ostringstream oss;
    oss << locus.object_id() << " "
        << length << " "
        << (mutation_rate.get() ? mutation_rate->object_id() : "null");
    return oss.str();
}


MutationGenerator_Regions::MutationInfos MutationGenerator_Regions::generate_mutations(
    const Population& population,
    size_t generation_index,
    size_t population_index) const
{
    MutationInfos result;

    for (RegionInfos::const_iterator region_info=region_infos_.begin(); region_info!=region_infos_.end(); ++region_info)
    {
        const double mu = region_info->mutation_rate->value(generation_index, population_index); // per site, per chromosome, per generation
        const size_t total_site_count = region_info->length * 2 * population.population_size();
        const double rate = mu * total_site_count;
        Random::DistributionPtr poisson = Random::create_poisson_distribution("id_dummy", rate);

        size_t mutant_count = size_t(poisson->random_value());

        vector<size_t> total_site_indices = 
            Random::random_indices_without_replacement(total_site_count, mutant_count);

        for (vector<size_t>::const_iterator index=total_site_indices.begin(); index!=total_site_indices.end(); ++index)
        {
            size_t position_index = *index % region_info->length;
            size_t chromosome_index = *index / region_info->length;

            MutationInfo info;
            info.individual_index = chromosome_index >> 1;
            info.locus.chromosome_pair_index = region_info->locus.chromosome_pair_index;
            info.locus.position = region_info->locus.position + position_index;
            info.which = chromosome_index & 0x01;
            info.value = 1;
            result.push_back(info);
        }
    }

    return result;
}


Parameters MutationGenerator_Regions::parameters() const
{
    Parameters parameters;
    for (RegionInfos::const_iterator it=region_infos_.begin(); it!=region_infos_.end(); ++it)
        parameters.insert_name_value("locus:length:rate", it->configuration());
    return parameters;
}


void MutationGenerator_Regions::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("locus:length:rate");

    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        istringstream iss(*it);
        string id_locus, id_trajectory;
        double length = 0;
        iss >> id_locus >> length >> id_trajectory;

        if (id_locus.empty() || id_trajectory.empty() || !iss)
            throw runtime_error(("[MutationGenerator_Regions] Error parsing locus:length:rate: " + *it).c_str());

        region_infos_.push_back(RegionInfo(*registry.get<Locus>(id_locus),
                                          static_cast<size_t>(length),
                                          registry.get<Trajectory>(id_trajectory)));
    }
}


void MutationGenerator_Regions::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (RegionInfos::const_iterator it=region_infos_.begin(); it!=region_infos_.end(); ++it)
    {
        it->locus.write_configuration(os, ids_written);
        it->mutation_rate->write_configuration(os, ids_written);
    }
}


