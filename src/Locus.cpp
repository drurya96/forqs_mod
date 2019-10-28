//
// Locus.cpp
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


#include "Locus.hpp"
#include "Random.hpp"
#include "Simulator.hpp"
#include <numeric>


using namespace std;


//
// Locus
//


Parameters Locus::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("chromosome", chromosome_pair_index + 1);
    parameters.insert_name_value("position", position);
    return parameters;
}


void Locus::configure(const Parameters& parameters, const Registry& registry)
{
    size_t chromosome = parameters.value<size_t>("chromosome", 1); // 1-based
    if (chromosome == 0) throw runtime_error("[Locus] Invalid chromosome == 0.");
    chromosome_pair_index = chromosome - 1;
    position = parameters.value<unsigned int>("position", 0);
}


bool operator<(const Locus& a, const Locus& b)
{
    return (a.chromosome_pair_index < b.chromosome_pair_index) ||
           (a.chromosome_pair_index == b.chromosome_pair_index && a.position < b.position);
}


bool operator==(const Locus& a, const Locus& b)
{
    return (a.chromosome_pair_index == b.chromosome_pair_index && a.position == b.position);
}


bool operator!=(const Locus& a, const Locus& b)
{
    return !(a==b);
}


std::ostream& operator<<(std::ostream& os, const Locus& locus)
{
    os << "(" << locus.chromosome_pair_index + 1 << "," << locus.position << ")";
    return os;
}


//
// LocusList
//


Parameters LocusList::parameters() const
{
    if (empty()) return Parameters();

    Parameters parameters;

    const_iterator it = begin();

    // unnamed loci

    for (; it!=end(); ++it)
    {
        size_t index = it - begin();

        if (it->object_id() == generate_locus_id(index))
        {
            ostringstream value;
            value << it->chromosome_pair_index + 1 << " " << it->position; // 1-based
            parameters.insert_name_value("chromosome:position", value.str());
        }
        else
        {
            break;
        }
    }

    // named loci

    ostringstream named_loci;
    for (; it!=end(); ++it)
        named_loci << it->object_id() << " ";
    if (!named_loci.str().empty())
        parameters.insert_name_value("loci", named_loci.str());

    return parameters;
}


void LocusList::configure(const Parameters& parameters, const Registry& registry)
{
    clear();

    vector<string> unnamed = parameters.values<string>("chromosome:position");
    for (vector<string>::const_iterator it=unnamed.begin(); it!=unnamed.end(); ++it)
    {
        size_t index = it - unnamed.begin();

        istringstream iss(*it);
        size_t chromosome = 0; // 1-based
        unsigned int position = 0;
        iss >> chromosome >> position;

        if (chromosome == 0)
            throw runtime_error("[LocusList] Invalid chromosome number.");

        push_back(Locus(generate_locus_id(index), chromosome-1, position));
    }

    vector<string> locus_ids = parameters.value_vector<string>("loci", vector<string>());
    for (vector<string>::const_iterator id=locus_ids.begin(); id!=locus_ids.end(); ++id)
    {
        LocusPtr locus = registry.get<Locus>(*id, nothrow);
        LocusListPtr locus_list = registry.get<LocusList>(*id, nothrow);

        if (locus.get()) 
            push_back(*locus);
        else if (locus_list.get())
            copy(locus_list->begin(), locus_list->end(), back_inserter(*this));
        else
            throw runtime_error("[LocusList] id must be Locus or LocusList.");
    }
}


void LocusList::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    // write only named loci

    for (const_iterator it=begin(); it!=end(); ++it)
    {
        size_t index = it - begin();
        if (it->object_id() == generate_locus_id(index)) continue;
        it->write_configuration(os, ids_written);
    }
}


string LocusList::generate_locus_id(size_t index) const
{
    ostringstream id;
    id << object_id() << "[" << index << "]";
    return id.str();
}


//
// LocusList_Random
//


Parameters LocusList_Random::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("locus_count", locus_count_);
    return parameters;
}


void LocusList_Random::configure(const Parameters& parameters, const Registry& registry)
{
    locus_count_ = parameters.value<size_t>("locus_count");
}


void LocusList_Random::initialize(const SimulatorConfig& config)
{
    // get chromosome lengths from PCG

    const size_t chromosome_pair_count = config.population_config_generator->chromosome_pair_count();
    const vector<unsigned int>& chromosome_lengths = config.population_config_generator->chromosome_lengths();

    if (chromosome_pair_count != chromosome_lengths.size())
    {
        ostringstream oss;
        oss << "[LocusList_Random] Number of chromosome_lengths " 
            << "doesn't match chromosome_pair_count in "
            << config.population_config_generator->class_name() << " "
            << config.population_config_generator->object_id();
        throw runtime_error(oss.str().c_str());
    }

    // generate random positions

    vector<long> cdf;
    long genome_length;

    cdf.resize(chromosome_lengths.size());
    partial_sum(chromosome_lengths.begin(), chromosome_lengths.end(), cdf.begin(), plus<long>());
    genome_length = cdf.empty() ? 0 : cdf.back();

    if (genome_length == 0)
        throw runtime_error("[LocusList_Random] Genome length 0.");

    for (size_t i=0; i<locus_count_; ++i)
    {
        long roll = long(Random::uniform_long(0, genome_length));

        size_t chromosome_index = lower_bound(cdf.begin(), cdf.end(), roll) - cdf.begin();

        size_t position = roll;
        if (chromosome_index > 0) position -= cdf[chromosome_index - 1];

        push_back(Locus(generate_locus_id(i), chromosome_index, position));
    }
}


