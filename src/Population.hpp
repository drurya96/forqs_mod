//
// Population.hpp
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


#ifndef _POPULATION_HPP_
#define _POPULATION_HPP_


#include "DataVector.hpp"
#include "PopulationData.hpp"
#include "ChromosomePairRange.hpp"
#include "shared_ptr.hpp"
#include <vector>


class MatingDistribution
{
    public:

    struct Entry
    {
        double weight;
        size_t first;
        std::string first_fitness;
        size_t second;
        std::string second_fitness;
        mutable bool valid;

        Entry(double _weight = 0, size_t _first = 0, size_t _second = 0)
        :   weight(_weight), first(_first), second(_second), valid(true)
        {}

        Entry(double _weight, 
              size_t _first, const std::string& _first_fitness, 
              size_t _second, const std::string& _second_fitness)
        :   weight(_weight), 
            first(_first), first_fitness(_first_fitness), 
            second(_second), second_fitness(_second_fitness), valid(true)
        {}
    };

    typedef std::vector<Entry> Entries;

    std::string default_fitness_function;
    void push_back(const Entry& entry);

    const Entries& entries() const {return entries_;} 
    const std::vector<double>& cumulative_weights() const {return cumulative_weights_;}
    bool empty() const {return entries_.empty();} 

    void validate_entries(const PopulationDataPtrs& population_datas) const;
    const Entry& random() const;

    private:

    Entries entries_;
    std::vector<double> cumulative_weights_;

    double total_weight() const {return cumulative_weights_.empty() ? 0 : cumulative_weights_.back();}
};


bool operator==(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b);
bool operator!=(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b);
bool operator==(const MatingDistribution& a, const MatingDistribution& b);
bool operator!=(const MatingDistribution& a, const MatingDistribution& b);
std::ostream& operator<<(std::ostream& os, const MatingDistribution::Entry& entry);
std::ostream& operator<<(std::ostream& os, const MatingDistribution& md);
std::istream& operator>>(std::istream& is, MatingDistribution::Entry& entry);
std::istream& operator>>(std::istream& is, MatingDistribution& md);


class Population;
typedef boost::shared_ptr<Population> PopulationPtr;
typedef std::vector<PopulationPtr> PopulationPtrs;
typedef boost::shared_ptr<PopulationPtrs> PopulationPtrsPtr;


class Population
{
    public:

    struct Config
    {
        size_t population_size;

        // for generating population from nothing
        size_t chromosome_pair_count; 
        unsigned int id_offset;

        // for generating population from a previous generation
        MatingDistribution mating_distribution;

        Config() 
        :   population_size(0), chromosome_pair_count(0), id_offset(0)
        {}
    };

    typedef std::vector<Config> Configs;

    size_t population_size() const {return population_size_;}
    size_t chromosome_pair_count() const {return chromosome_pair_count_;}
    bool empty() const {return (population_size_ == 0);}

    void read_text(std::istream& is);
    void write_text(std::ostream& os) const;
    void read_binary(std::istream& is);
    void write_binary(std::ostream& os) const;

    void create_organisms(const Config& config,
                          const PopulationPtrs& populations,
                          const PopulationDataPtrs& population_datas,
                          const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array);

    // convenience function: creates new generation from previous by calling create_organisms() for each Population

    static PopulationPtrsPtr create_populations(const Configs& configs,
                                                const PopulationPtrs& previous, 
                                                const PopulationDataPtrs& population_datas, 
                                                const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array);

    // implementation-dependent range iteration

    virtual void allocate_memory() = 0;

    virtual ChromosomePairRangeIterator begin() = 0;
    virtual const ChromosomePairRangeIterator begin() const = 0;

    virtual ChromosomePairRangeIterator end() = 0;
    virtual const ChromosomePairRangeIterator end() const = 0;

    virtual ChromosomePairRange chromosome_pair_range(size_t organism_index) = 0;
    virtual const ChromosomePairRange chromosome_pair_range(size_t organism_index) const = 0;

    virtual ~Population() {}

    protected:

    Population()
    :   population_size_(0), chromosome_pair_count_(0)
    {}

    size_t population_size_;
    size_t chromosome_pair_count_;

    private:

    // disallow copying
    Population(Population&);
    Population& operator=(Population&);
};


bool operator==(const Population::Config& a, const Population::Config& b);
bool operator!=(const Population::Config& a, const Population::Config& b);
std::ostream& operator<<(std::ostream& os, const Population::Config& config);
std::istream& operator>>(std::istream& is, Population::Config& config);

// multiple generations: each generation is Population::Configs, one Config per population
std::ostream& operator<<(std::ostream& os, const std::vector<Population::Configs>& generations);
std::istream& operator>>(std::istream& is, std::vector<Population::Configs>& generations);


bool operator==(const Population& a, const Population& b);
bool operator!=(const Population& a, const Population& b);
std::ostream& operator<<(std::ostream& os, const Population& p);    // calls write_text()
std::istream& operator>>(std::istream& is, Population& p);          // calls read_text()


#endif // _POPULATION_HPP_

