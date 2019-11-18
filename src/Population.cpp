//
// Population.cpp
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


#include "Population.hpp"
#include "Population_ChromosomePairs.hpp"
#include "Random.hpp"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iterator>
#include <set>
#include <fstream>
#include <cmath>


using namespace std;


//
// MatingDistribution
//


void MatingDistribution::push_back(const Entry& entry)
{
	entries_.push_back(entry);
	cumulative_weights_.push_back(total_weight() + entry.weight);
}


void MatingDistribution::validate_entries(const PopulationDataPtrs& population_datas) const
{
	bool something_valid = false;

	for (Entries::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
	{
		const string& ff1 = !it->first_fitness.empty() ? it->first_fitness : default_fitness_function;
		const string& ff2 = !it->second_fitness.empty() ? it->second_fitness : default_fitness_function;
		it->valid = ((ff1.empty() || 
						!population_datas.at(it->first)->trait_values->get(ff1)->all_zero()) && 
					 (ff2.empty() || 
						!population_datas.at(it->second)->trait_values->get(ff2)->all_zero()));
		if (it->valid) something_valid = true;

		if (!it->valid)
			cerr << "[MatingDistribution] Warning: marking invalid entry " << *it << endl;
	}

	if (!something_valid)
		throw runtime_error("[MatingDistribution] No valid entries.");
}


const MatingDistribution::Entry& MatingDistribution::random() const
{
	const size_t max_attempts = 10000;

	for (size_t attempt=0; attempt<max_attempts; ++attempt)
	{
		double roll = Random::uniform_real(0, total_weight());

		vector<double>::const_iterator it = 
			lower_bound(cumulative_weights_.begin(), cumulative_weights_.end(), roll);

		if (it == cumulative_weights_.end())
		{
			cout << "total weight: " << total_weight() << endl;
			cout << "roll: " << roll << endl;
			throw runtime_error("[MatingDistribution::random()] This isn't happening.");
		}

		const Entry& entry = entries_[it - cumulative_weights_.begin()];
		
		if (entry.valid) return entry;
	}

	throw runtime_error("[MatingDistribution] Something's wrong: no valid entries?");
}


bool operator==(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b)
{
	return (a.first == b.first && 
			a.first_fitness == b.first_fitness && 
			a.second == b.second &&
			a.second_fitness == b.second_fitness);
}


bool operator!=(const MatingDistribution::Entry& a, const MatingDistribution::Entry& b)
{
	return !(a==b);
}


bool operator==(const MatingDistribution& a, const MatingDistribution& b)
{
	if (a.entries().size() != b.entries().size()) return false;

	const double epsilon = 1e-12;

	for (size_t i=0; i<a.entries().size(); ++i)
		if (a.entries()[i] != b.entries()[i] ||
			fabs(a.cumulative_weights()[i] - b.cumulative_weights()[i]) > epsilon)
			return false;

	return true;
}


bool operator!=(const MatingDistribution& a, const MatingDistribution& b)
{
	return !(a==b);
}


ostream& operator<<(ostream& os, const MatingDistribution::Entry& entry)
{
	// <.6|0,0>
	os << "<" << entry.weight << "|" 
	   << entry.first << "," 
	   << entry.first_fitness << "," 
	   << entry.second  << ","
	   << entry.second_fitness << ">";
	return os;
}


ostream& operator<<(ostream& os, const MatingDistribution& md)
{
	// {<.6|0,0><.2|1,0><.1|2,0>}
	os << "{";
	copy(md.entries().begin(), md.entries().end(), ostream_iterator<MatingDistribution::Entry>(os,""));
	os << "}";
	return os;
}


istream& operator>>(istream& is, MatingDistribution::Entry& entry)
{
	// <0.2|1,ff1,11,ff2>

	string buffer;
	getline(is, buffer, '>');
	if (!is) return is;

	char open, pipe, comma1, comma3;

	istringstream iss(buffer);
	iss >> open >> entry.weight >> pipe;
	iss >> entry.first >> comma1;
	getline(iss, entry.first_fitness, ','); // comma2
	iss >> entry.second >> comma3 >> entry.second_fitness;

	if (!is || open!='<' || pipe!='|' || comma1!=',' || comma3!=',')
		throw runtime_error("[operator>>(MatingDistribution::Entry)] Invalid input string.");

	return is;
}


istream& operator>>(istream& is, MatingDistribution& md)
{
	string buffer;
	getline(is, buffer,'}');
	if (!is) return is;

	char open;
	istringstream iss(buffer);
	iss >> open;
	if (open!='{')
		throw runtime_error("[operator>>(MatingDistribution)] Invalid input string.");

	// read in vector of entries

	MatingDistribution::Entries entries;
	copy(istream_iterator<MatingDistribution::Entry>(iss), 
		 istream_iterator<MatingDistribution::Entry>(),
		 back_inserter(entries));	 

	// clear any old data

	md = MatingDistribution();

	// push_back() does the cumulative weight calculation

	for (MatingDistribution::Entries::const_iterator it=entries.begin(); it!=entries.end(); ++it)
		md.push_back(*it);

	return is;
}


//
// Population::Config
//


bool operator==(const Population::Config& a, const Population::Config& b)
{
	return a.population_size == b.population_size &&
		   a.id_offset == b.id_offset &&
		   a.chromosome_pair_count == b.chromosome_pair_count &&
		   a.mating_distribution == b.mating_distribution;
}


bool operator!=(const Population::Config& a, const Population::Config& b)
{
	return !(a==b);
}


ostream& operator<<(ostream& os, const Population::Config& config)
{
	os << "population"
	   << " population_size=" << config.population_size;

	if (config.id_offset != 0)
		os << " id_offset=" << config.id_offset;

	if (config.chromosome_pair_count != 0)
		os << " chromosome_pair_count=" << config.chromosome_pair_count;

	if (!config.mating_distribution.entries().empty())
		os << " mating_distribution=" << config.mating_distribution;

	return os;
}


istream& operator>>(istream& is, Population::Config& config)
{
	// population population_size=4 id_offset=1000 chromosome_pair_count=3 mating_distribution={<0.42|0,0><0.66|1,0><0.23|2,1>}

	string buffer;
	getline(is, buffer);

	vector<string> tokens;
	istringstream iss(buffer);
	copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

	if (tokens.empty() || tokens[0]!="population")
	{
		cerr << "--> " << buffer << endl;
		throw runtime_error("[operator<<(Population::Config)] Invalid config format");
	}

	for (vector<string>::const_iterator it=tokens.begin()+1; it!=tokens.end(); ++it)
	{
		size_t index_equal = it->find('=');
		if (index_equal == string::npos)
		{
			cerr << "Ignoring invalid token: " << *it << endl;
			continue;
		}

		string name = it->substr(0,index_equal);
		istringstream value(it->substr(index_equal+1));

		if (name.empty() || value.str().empty())
		{
			cerr << "Ignoring invalid token: " << *it << endl;
			continue;
		}

		if (name == "population_size")
			value >> config.population_size;
		else if (name == "id_offset")
			value >> config.id_offset;
		else if (name == "chromosome_pair_count")
			value >> config.chromosome_pair_count;
		else if (name == "mating_distribution")
			value >> config.mating_distribution;
	}

	return is;
}


ostream& operator<<(ostream& os, const vector<Population::Configs>& generations)
{
	for (size_t i=0; i<generations.size(); i++)
	{
		os << "generation " << i << endl;
		copy(generations[i].begin(), generations[i].end(), 
			 ostream_iterator<Population::Config>(os,"\n"));
		os << endl;
	}

	return os;
}


istream& operator>>(istream& is, vector<Population::Configs>& generations)
{
	size_t current_id_offset = 0;

	while (is)
	{
		// parse line by line

		string buffer;
		getline(is, buffer);
		if (!is) return is;

		vector<string> tokens;
		istringstream iss(buffer);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

		// switch on first token

		if (tokens.empty() || !tokens[0].empty() && tokens[0][0]=='#')
			continue;
		else if (tokens[0] == "generation")
		{
			generations.push_back(vector<Population::Config>());
			current_id_offset = 0;
		}
		else if (tokens[0] == "population")
		{
			generations.back().push_back(Population::Config());
			istringstream temp(buffer);
			Population::Config& config = generations.back().back();
			temp >> config;

			// automatically set id_offset if not specified
			if (config.chromosome_pair_count > 0) // initial generation
			{
				if (config.id_offset == 0)
					config.id_offset = current_id_offset;
				current_id_offset += 2*config.population_size;
			}
		}
		else
			cerr << "Ignoring invalid configuration line:\n" << buffer << endl;
	}

	return is;
}


//
// Population
//


void Population::read_text(std::istream& is)
{
	// read header lines

	string buffer;

	getline(is, buffer);
	istringstream iss1(buffer);
	string population_size_string;
	iss1 >> population_size_string >> population_size_;

	getline(is, buffer);
	istringstream iss2(buffer);
	string chromosome_pair_count_string;
	iss2 >> chromosome_pair_count_string >> chromosome_pair_count_;

	getline(is, buffer);

	if (population_size_string != "population_size" ||
		chromosome_pair_count_string != "chromosome_pair_count")
		throw runtime_error("[Population::read_text] Bad Population text format.");

	// read organisms

	allocate_memory();

	for (ChromosomePairRangeIterator it=begin(); it!=end(); ++it)
	{
		for (ChromosomePair* p=it->begin(); p!=it->end(); ++p)
		{
			getline(is, buffer);
			istringstream iss_plus(buffer);
			char plus;
			iss_plus >> plus >> p->first;
			if (plus != '+') throw runtime_error("[Population::read_text()] + expected");

			getline(is, buffer);
			istringstream iss_minus(buffer);
			char minus;
			iss_minus >> minus >> p->second;
			if (minus != '-') throw runtime_error("[Population::read_text()] - expected");
		}
		
		getline(is, buffer);
	}
}


void Population::write_text(std::ostream& os) const
{
	// write header lines

	os << "population_size " << population_size_ << endl;
	os << "chromosome_pair_count " << chromosome_pair_count_ << endl;
	os << endl;

	// write organisms

	for (const ChromosomePairRangeIterator it=begin(); it!=end(); ++it)
	{
		for (const ChromosomePair* p=it->begin(); p!=it->end(); ++p)
			os << "+ " << p->first << "\n- " << p->second << endl;
		os << endl;
	}
}


void Population::read_binary(istream& is)
{
	// read header

	is.read((char*)&population_size_, sizeof(size_t));
	is.read((char*)&chromosome_pair_count_, sizeof(size_t));

	// read organisms

	allocate_memory();

	for (ChromosomePairRangeIterator it=begin(); it!=end(); ++it)
	{
		for (ChromosomePair* p=it->begin(); p!=it->end(); ++p)
		{
			p->first.read(is);
			p->second.read(is);
		}
	}
}


void Population::write_binary(ostream& os) const
{
	// write header

	os.write((const char*)&population_size_, sizeof(size_t));
	os.write((const char*)&chromosome_pair_count_, sizeof(size_t));

	// write organisms

	for (ChromosomePairRangeIterator it=begin(); it!=end(); ++it)
	{
		for (ChromosomePair* p=it->begin(); p!=it->end(); ++p)
		{
			p->first.write(os);
			p->second.write(os);
		}
	}
}


namespace {


class RandomOrganismIndexGenerator
{
	public:

	RandomOrganismIndexGenerator(size_t population_size,
								 const DataVectorPtr& fitness_vector)
	:   population_size_(population_size),
		fitness_cdf_max_(0)
	{
		if (population_size == 0)
			throw runtime_error("[RandomOrganismIndexGenerator] Population size 0.");

		if (fitness_vector.get()) 
		{
			fitness_cdf_ = fitness_vector->cdf(); // memory allocation for cdf

			if (!fitness_cdf_.get() || fitness_cdf_->empty() || fitness_cdf_->size() != population_size_)
				throw runtime_error("[RandomOrganismIndexGenerator] This isn't happening.");

			fitness_cdf_max_ = fitness_cdf_->back();
		}
	}

	size_t operator()() const
	{
		if (!fitness_cdf_.get()) 
		{
			return Random::uniform_integer(0, population_size_-1); // uniform random index
		}
		else
		{
			// pick random index according to fitnesses
			double roll = Random::uniform_real(0, fitness_cdf_max_);
			DataVector::const_iterator it = lower_bound(fitness_cdf_->begin(), fitness_cdf_->end(), roll);
			return it - fitness_cdf_->begin();
		}
	}

	private:

	size_t population_size_;
	DataVectorPtr fitness_cdf_;
	double fitness_cdf_max_;

	//allows Population (need Population::create_populations) to access private members (A.D. 11/2019)
	friend void Population::create_organisms(const Config& config,
								  const PopulationPtrs& populations,
								  const PopulationDataPtrs& population_datas,
								  const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array);
	
};


typedef boost::shared_ptr<RandomOrganismIndexGenerator> RandomOrganismIndexGeneratorPtr;


class RandomOrganismIndexGeneratorMap
{
	public:

	RandomOrganismIndexGeneratorMap(const PopulationDataPtrs& population_datas,
									string default_fitness_function)
	:   population_datas_(population_datas),
		default_fitness_function_(default_fitness_function)
	{}

	RandomOrganismIndexGeneratorPtr get(size_t population_index, string fitness_function)
	{
		if (fitness_function.empty()) 
			fitness_function = default_fitness_function_;
		
		Key key = make_pair(population_index, fitness_function);

		if (!generator_map_.count(key)) // lazy initialization of generator
		{
			if (population_index >= population_datas_.size())
				throw runtime_error("[Population::RandomOrganismIndexGeneratorMap] Bad index.");

			const PopulationData& popdata = *population_datas_[population_index];
			size_t population_size = popdata.population_size;

			DataVectorPtr fitness = fitness_function.empty() ?  DataVectorPtr() : 
				popdata.trait_values->get(fitness_function);

			if (fitness.get() && fitness->size() != population_size)
				throw runtime_error("[Population::RandomOrganismIndexGeneratorMap] Bad population size.");

			generator_map_[key] =
				RandomOrganismIndexGeneratorPtr(
					new RandomOrganismIndexGenerator(population_size, fitness));
		}

		return generator_map_[key];
	}

	private:

	const PopulationDataPtrs& population_datas_;
	string default_fitness_function_;

	typedef pair<size_t,string> Key;
	typedef map<Key,RandomOrganismIndexGeneratorPtr> GeneratorMap;
	GeneratorMap generator_map_;
};


} // namespace


void Population::create_organisms(const Config& config,
								  const PopulationPtrs& populations,
								  const PopulationDataPtrs& population_datas,
								  const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array)
{
	if (config.population_size == 0)
		return;

	// sanity check: chromosome_pair_count must be specified

	if (config.chromosome_pair_count==0)
		throw runtime_error("[Population::create_organisms()] Must specify nonzero chromosome pair count.");

	// allocate memory

	population_size_ = config.population_size;
	chromosome_pair_count_ = config.chromosome_pair_count;
	allocate_memory();

	// create organisms from nothing

	if (populations.empty())
	{
		if (config.chromosome_pair_count == 0)
			throw runtime_error("[Population::create_organisms()] Chromosome pair count 0.\n");

		ChromosomePairRangeIterator range = begin();
		
		for (size_t i=0; i<config.population_size; ++i, ++range)
		{
			unsigned int id0 = config.id_offset + 2*i;
			range->create_child(id0, id0+1); 
		}

		return;
	}

	// create organisms from previous generation

	// sanity checks

	if (populations.size() != population_datas.size())
		throw runtime_error("[Population::create_organisms()] Population count and data count differ.");

	PopulationDataPtrs::const_iterator popdata = population_datas.begin();
	for (PopulationPtrs::const_iterator population = populations.begin();
		 population != populations.end(); ++population, ++popdata)
	{
		if ((*population)->population_size() != (*popdata)->population_size)
			throw runtime_error("[Population::create_organisms()] Population size mismatch.");
	}

	for (MatingDistribution::Entries::const_iterator it=config.mating_distribution.entries().begin();
		 it!=config.mating_distribution.entries().end(); ++it)
	{
		if (max(it->first, it->second) >= populations.size())
			throw runtime_error("[Population::create_organisms()] Indices out of bounds.");
	}

	config.mating_distribution.validate_entries(population_datas);

	// instantiate RandomOrganismIndexGeneratorMap

	RandomOrganismIndexGeneratorMap generator_map(population_datas,
		config.mating_distribution.default_fitness_function);

	// create Organisms for new population

	ChromosomePairRangeIterator range_child = begin();

	for (size_t i=0; i<config.population_size; ++i, ++range_child)
	{

		const MatingDistribution::Entry& entry = config.mating_distribution.random();

		const RandomOrganismIndexGenerator& generator_mom = 
			*generator_map.get(entry.first, entry.first_fitness);

		const RandomOrganismIndexGenerator& generator_dad = 
			*generator_map.get(entry.second, entry.second_fitness);

		size_t index_mom = generator_mom();
		size_t index_dad = 0;
		do { // avoid selfing
			index_dad = generator_dad();
		} while (entry.first == entry.second && index_mom == index_dad);

		ChromosomePairRange range_mom = populations[entry.first]->chromosome_pair_range(index_mom);
		ChromosomePairRange range_dad = populations[entry.second]->chromosome_pair_range(index_dad);
		

		if (recombination_position_generators_array.size() != 1){
			// These values are the trait values of the parents chosen to mate
			double mom_value = generator_mom.fitness_cdf_->at(index_mom);
			double dad_value = generator_dad.fitness_cdf_->at(index_dad);

			range_mom.recombination_rate = mom_value*100000;
			range_dad.recombination_rate = dad_value*100000;
		}

		range_child->create_child(range_mom, range_dad, recombination_position_generators_array);
	}
}


// static
PopulationPtrsPtr Population::create_populations(const Population::Configs& configs,
												 const PopulationPtrs& previous, 
												 const PopulationDataPtrs& population_datas,
												 const RecombinationPositionGeneratorPtrsArray& recombination_position_generators_array)
{
	PopulationPtrsPtr result(new PopulationPtrs);

	for (vector<Population::Config>::const_iterator it=configs.begin(); it!=configs.end(); ++it)
	{
		PopulationPtr p(new Population_ChromosomePairs);
		p->create_organisms(*it, previous, population_datas, recombination_position_generators_array);
		result->push_back(p);
	}		

	return result;
}


//
// Population operators
//


bool operator==(const Population& a, const Population& b)
{
	if (a.population_size() != b.population_size()) return false;
	if (a.chromosome_pair_count() != b.chromosome_pair_count()) return false;

	const ChromosomePairRangeIterator range_a_end = a.end();
	for (const ChromosomePairRangeIterator range_a=a.begin(), range_b=b.begin();
		 range_a != range_a_end; ++range_a, ++range_b)
	{
		if (range_a->size() != range_b->size()) return false;

		for (const ChromosomePair* p_a=range_a->begin(), * p_b=range_b->begin(); 
			 p_a!=range_a->end(); ++p_a, ++p_b)
		{
			if (*p_a != *p_b)
			{
				//cout << "bad: " << p_a->first << " " << p_a->second << " " << p_b->first << " " << p_b->second << endl;
				return false;
			}
		}
	}

	return true;
}


bool operator!=(const Population& a, const Population& b)
{
	return !(a==b);
}


ostream& operator<<(ostream& os, const Population& p)
{
	p.write_text(os);
	return os;
}


istream& operator>>(istream& is, Population& p)
{
	p.read_text(is);
	return is;
}


