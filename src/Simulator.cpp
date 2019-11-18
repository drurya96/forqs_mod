//
// Simulator.cpp
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


#include "Simulator.hpp"
#include <iostream>
#include <iterator>
#include <cmath>
#include <algorithm>

#include "Population_ChromosomePairs.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"



using namespace std;


//
// SimulatorConfig
//


SimulatorConfig::SimulatorConfig(const string& id)
:   Configurable(id), 
	seed(0), 
	write_popconfig(false),
	write_vi(false),
	use_random_seed(false)
{}


Parameters SimulatorConfig::parameters() const
{
	Parameters parameters;

	parameters.insert_name_value("seed", seed);
	parameters.insert_name_value("output_directory", output_directory);
	parameters.insert_name_value("write_popconfig", write_popconfig);
	parameters.insert_name_value("write_vi", write_vi);

	if (population_config_generator.get())
		parameters.insert_name_value("population_config_generator", population_config_generator->object_id());

	for (RecombinationPositionGeneratorPtrs::const_iterator it=recombination_position_generators.begin();
		it!=recombination_position_generators.end(); ++it)
		if (it->get())
			parameters.insert_name_value("recombination_position_generator", (*it)->object_id());

	//make_recombination_generators(parameters);

	if (variant_indicator.get())
		parameters.insert_name_value("variant_indicator", variant_indicator->object_id());

	for (vector<QuantitativeTraitPtr>::const_iterator it=quantitative_traits.begin(); it!=quantitative_traits.end(); ++it)
		if (it->get())
			parameters.insert_name_value("quantitative_trait", (*it)->object_id());

	if (mutation_generator.get())
		parameters.insert_name_value("mutation_generator", mutation_generator->object_id());

	for (vector<ReporterPtr>::const_iterator it=reporters.begin(); it!=reporters.end(); ++it)
		if (it->get())
			parameters.insert_name_value("reporter", (*it)->object_id());

	return parameters;
}


void SimulatorConfig::configure(const Parameters& parameters, const Registry& registry)
{
	if (parameters.count("seed"))
		seed = parameters.value<unsigned int>("seed");
	else
		use_random_seed = true;

	this->registry = registry;
	output_directory = parameters.value<string>("output_directory", "");
	write_popconfig = parameters.value<bool>("write_popconfig", false);
	write_vi = parameters.value<bool>("write_vi", false);

	population_config_generator = registry.get<PopulationConfigGenerator>(
		parameters.value<string>("population_config_generator"));

	vector<string> rpg_ids = parameters.values<string>("recombination_position_generator");
	for (vector<string>::const_iterator it=rpg_ids.begin(); it!=rpg_ids.end(); ++it){
		recombination_position_generators.push_back(
			registry.get<RecombinationPositionGenerator>(*it));
	}

	string mutation_generator_id = parameters.value<string>("mutation_generator", "");
	if (!mutation_generator_id.empty())
		mutation_generator = registry.get<MutationGenerator>(mutation_generator_id);

	string variant_indicator_id = parameters.value<string>("variant_indicator", "");
	if (!variant_indicator_id.empty())
		variant_indicator = registry.get<VariantIndicator>(variant_indicator_id);

	vector<string> qtids = parameters.values<string>("quantitative_trait");
	for (vector<string>::const_iterator it=qtids.begin(); it!=qtids.end(); ++it)
		quantitative_traits.push_back(registry.get<QuantitativeTrait>(*it));

	vector<string> reporter_ids = parameters.values<string>("reporter");
	for (vector<string>::const_iterator it=reporter_ids.begin(); it!=reporter_ids.end(); ++it){
		reporters.push_back(registry.get<Reporter>(*it));
	}
}

void SimulatorConfig::make_recombination_generators(Parameters& parameters) const
{

	//recombination_position_generators.clear();
	for (RecombinationPositionGeneratorPtrs::const_iterator it=recombination_position_generators.begin();
		it!=recombination_position_generators.end(); ++it)
		if (it->get())
			parameters.insert_name_value("recombination_position_generator", (*it)->object_id());


}

void SimulatorConfig::write_child_configurations(ostream& os, set<string>& ids_written) const
{
	if (population_config_generator.get())
		population_config_generator->write_configuration(os, ids_written);

	for (RecombinationPositionGeneratorPtrs::const_iterator it=recombination_position_generators.begin(); 
		it!=recombination_position_generators.end(); ++it)
		if (it->get())
			(*it)->write_configuration(os, ids_written);

	if (variant_indicator.get())
		variant_indicator->write_configuration(os, ids_written);

	for (vector<QuantitativeTraitPtr>::const_iterator it=quantitative_traits.begin(); it!=quantitative_traits.end(); ++it)
		if (it->get())
			(*it)->write_configuration(os, ids_written);

	if (mutation_generator.get())
		mutation_generator->write_configuration(os, ids_written);

	for (vector<ReporterPtr>::const_iterator it=reporters.begin(); it!=reporters.end(); ++it)
		if (it->get())
			(*it)->write_configuration(os, ids_written);
}


//
// Simulator
//


Simulator::Simulator(SimulatorConfig& config, const Parameters& parameters)
:   config_(config),
	current_generation_index_(0), 
	current_populations_(new PopulationPtrs),
	current_population_datas_(new PopulationDataPtrs),
	update_step_(1)
{
	//cout << "1st in Simulator: config_.recombination_position_generators size: " << config_.recombination_position_generators.size() << endl;
	const size_t generation_count = config_.population_config_generator->generation_count();
	update_step_ = max(int(pow(10.0, int(log10(generation_count))-1)), 1);
	if (update_step_ > 10000) update_step_ = 10000;

	//A.D. 11/07/2019

	config_.recombination_position_generators_array.resize(1);
	config_.recombination_position_generators_array[0] = config_.recombination_position_generators;

	string cname;
	//Parameters temp_params;
	for (Configurable::Registry::iterator it = config_.registry.begin(); it != config_.registry.end(); it++){
		cname = (it->second)->class_name();
		if(cname == "FitnessFunction_Recombination"){
			FitnessFunction_Recombination * ffptr = dynamic_cast<FitnessFunction_Recombination*>(&(*(it->second)));
			FitnessFunction_Recombination ff = *ffptr;
			cout << ff.modify_status() << endl;
			if (ff.modify_status()){
				set_recombination_generators(ff, config_.recombination_position_generators_array);
			}
		}
	}
}


void set_recombination_generators(FitnessFunction_Recombination& ff, RecombinationPositionGeneratorPtrsArray& rpg_vector){
	int number_of_indecies = 1000000;
	rpg_vector.resize(number_of_indecies);

	//RecombinationPositionGeneratorPtrs default_rpgptrs = rpg_vector[0];
	
	/*
	RecombinationPositionGeneratorPtrs default_rpgptrs;
	default_rpgptrs.resize(rpg_vector[0].size());
	
	for (unsigned int i=0; i < rpg_vector[0].size(); i++){
		// rpg_vector[0][i] is an rpg_pointer
		// *rpg_vector[0][i] is an rpg
		RecombinationPositionGenerator_Uniform* oldptr = dynamic_cast<RecombinationPositionGenerator_Uniform*>(&(*rpg_vector[0][i]));
		RecombinationPositionGenerator_Uniform old = *oldptr;
		//RecombinationPositionGenerator_Uniform old_copy = RecombinationPositionGenerator_Uniform(old);
		//RecombinationPositionGeneratorPtr old_copy_ptr = RecombinationPositionGeneratorPtr(old_copy);
		//default_rpgptrs[i] = old_copy_ptr;

		RecombinationPositionGeneratorPtr old_copy_ptr = RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Uniform(old));
		default_rpgptrs[i] = old_copy_ptr;
	}
	*/
	


	
	for(int i=0; i < number_of_indecies; i++){
		//rpgptrs is a vector of rpg pointers
		double new_rate = (double)i/number_of_indecies;
		(void) new_rate;
		//RecombinationPositionGeneratorPtrs new_rpgptrs = default_rpgptrs;

		RecombinationPositionGeneratorPtrs new_rpgptrs;
		new_rpgptrs.resize(rpg_vector[0].size());
		for (unsigned int i=0; i < rpg_vector[0].size(); i++){
			// rpg_vector[0][i] is an rpg_pointer
			// *rpg_vector[0][i] is an rpg
			RecombinationPositionGenerator_Uniform* oldptr = dynamic_cast<RecombinationPositionGenerator_Uniform*>(&(*rpg_vector[0][i]));
			RecombinationPositionGenerator_Uniform old = *oldptr;
			RecombinationPositionGeneratorPtr old_copy_ptr = RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Uniform(old));
			new_rpgptrs[i] = old_copy_ptr;
		}

		//step through new_rpgptrs
		RecombinationPositionGeneratorPtrs::iterator it;
		for (it = new_rpgptrs.begin(); it < new_rpgptrs.end(); it++){

			RecombinationPositionGenerator_Uniform* old_generatorptr = dynamic_cast<RecombinationPositionGenerator_Uniform*>(&(**it));
			//cout << typeid(old_generatorptr).name() << endl;
			RecombinationPositionGenerator_Uniform old_generator = *old_generatorptr;
			// it is a pointer to a rpg pointer
			// **it is a rpg

			// we need to step through chromsomeInfos......
			RecombinationPositionGenerator_Uniform::ChromosomeInfos new_chromsome_infos;

			RecombinationPositionGenerator_Uniform::ChromosomeInfos::iterator ft;

			for (ft = old_generator.infos_.begin(); ft != old_generator.infos_.end(); ft++){
				// Make a new ChromosomeInfo and add it to new Infos
				RecombinationPositionGenerator_Uniform::ChromosomeInfo new_chromosome_info((*ft).length, new_rate);
				new_chromsome_infos.push_back(new_chromosome_info);
			}

			// create a new rpg
			RecombinationPositionGenerator_Uniform new_rpg(old_generator.object_id(), new_chromsome_infos);
			(**it) = new_rpg;
		}

		if (i % 100000 == 0){cout << i << endl;}

		// add a pointer to the new rpg to the rpg vector
		rpg_vector[i] = new_rpgptrs;
		//rpg_vector[i] = rpg_vector[0];
	}

	cout << "lower bound is " << ff.getLowerBound() << endl;
	cout << "upper bound is " << ff.getUpperBound() << endl;
}


namespace {

void generate_mutations(const PopulationPtrs& next_populations, 
						const MutationGenerator& mutation_generator,
						VariantIndicator& variant_indicator,
						size_t current_generation_index)
{
	size_t next_population_count = next_populations.size();

	PopulationPtrs::const_iterator population = next_populations.begin();
	for (size_t population_index=0; population_index!=next_population_count; ++population_index, ++population)
	{
		MutationGenerator::MutationInfos mutation_infos = 
			mutation_generator.generate_mutations(**population, 
												  current_generation_index,
												  population_index);
	
		for (MutationGenerator::MutationInfos::const_iterator info=mutation_infos.begin(); info!=mutation_infos.end(); ++info)
		{
			ChromosomePair* chromosome_pair = (*population)->chromosome_pair_range(info->individual_index).begin()
				+ info->locus.chromosome_pair_index;

			Chromosome& chromosome = info->which ? chromosome_pair->second : chromosome_pair->first;
			
			HaplotypeChunk& chunk = *chromosome.find_haplotype_chunk(info->locus.position);
			
			unsigned int new_id = variant_indicator.mutate(chunk.id, info->locus, info->value);

			chunk.id = new_id;
		}
	}
}


Loci construct_loci_list(const QuantitativeTraitPtrs& quantitative_traits,
						 const ReporterPtrs& reporters,
						 size_t current_generation_index,
						 bool is_final_generation)
{
	Loci loci_all;

	for (QuantitativeTraitPtrs::const_iterator qt=quantitative_traits.begin();
		 qt!=quantitative_traits.end(); ++qt)
	{
		const Loci& loci = (*qt)->loci();
		for (Loci::const_iterator locus=loci.begin(); locus!=loci.end(); ++locus)
			loci_all.insert(*locus);
	}

	for (ReporterPtrs::const_iterator reporter=reporters.begin();
		 reporter!=reporters.end(); ++reporter)
	{
		Loci loci = (*reporter)->loci(current_generation_index, is_final_generation);
		for (Loci::const_iterator locus=loci.begin(); locus!=loci.end(); ++locus)
			loci_all.insert(*locus);
	}

	return loci_all;
}

} // namespace


void Simulator::simulate_single_generation() // main loop iteration
{
	// sanity checks

	const size_t generation_count = config_.population_config_generator->generation_count();

	if (!current_populations_.get())
		throw runtime_error("[Simulator::simulate_single_generation()] Null pointer.");

	if (current_generation_index_ > generation_count)
		throw runtime_error("[Simulator::simulate_single_generation()] Population config not specified.");

	if (current_generation_index_ % update_step_ == 0)
		cout << "[Simulator] Generation " << current_generation_index_ << endl;

	// create next generation

	Population::Configs popconfigs = 
		config_.population_config_generator->population_configs(current_generation_index_,*current_population_datas_);

	PopulationPtrsPtr next_populations = Population::create_populations(
		popconfigs, 
		*current_populations_, 
		*current_population_datas_, 
		config_.recombination_position_generators_array);

	// generate mutations

	if (config_.mutation_generator.get())
		generate_mutations(*next_populations,
						   *config_.mutation_generator, 
						   *config_.variant_indicator, 
						   current_generation_index_);

	// construct loci list

	const bool is_final_generation = (current_generation_index_ == generation_count);
	
	Loci loci_all = construct_loci_list(config_.quantitative_traits,
										config_.reporters,
										current_generation_index_,
										is_final_generation);

	// allocate PopulationData for each population
	// note: construct individually, since each one allocates memory for maps

	size_t next_population_count = next_populations->size();

	PopulationDataPtrsPtr next_population_datas(new PopulationDataPtrs);

	for (size_t i=0; i<next_population_count; ++i)
		next_population_datas->push_back(PopulationDataPtr(new PopulationData));

	// fill in population metadata and calculate genotypes

	PopulationPtrs::const_iterator population = next_populations->begin();
	PopulationDataPtrs::iterator popdata = next_population_datas->begin();
	PopulationDataPtrs::iterator currentpopdata = current_population_datas_->begin();
	for (size_t population_index=0; population_index!=next_population_count; 
		 ++population_index, ++population, ++popdata, ++currentpopdata)
	{

		// Somewhere here I need to try to get the ChromosomePairRanges into the PopulationData..............

		// Actually instead try inside Population::create_populations...................

		(*popdata)->generation_index = current_generation_index_;
		(*popdata)->population_index = population_index;
		(*popdata)->population_size = (*population)->population_size();
		//(*popdata)->organisms = (*currentpopdata)->organisms;

		genotyper_.genotype(loci_all, **population, *config_.variant_indicator, 
			*(*popdata)->genotypes);
	}

	// calculate quantitative trait values

	for (QuantitativeTraitPtrs::const_iterator qt=config_.quantitative_traits.begin();
		 qt!=config_.quantitative_traits.end(); ++qt)
	{
		(*qt)->calculate_trait_values(*next_population_datas);
	}

	// update reporters

	current_populations_ = next_populations;
	current_population_datas_ = next_population_datas;

	/*
	// ------------------------------------------------------------------------------------------------------------------
	// This was all just to check how next_populations works and is in no way functional --A.D.

	PopulationPtrs::iterator vecpopptrit = (*next_populations).begin();

	for (; vecpopptrit != (*next_populations).end(); vecpopptrit++){
		Population_ChromosomePairs * tmpPairPtr = dynamic_cast<Population_ChromosomePairs*>(&(**vecpopptrit));
		//cout << typeid(tmpPairPtr).name() << endl;
	}


	//cout << "ID of current_populations_: " << typeid(current_populations_).name() << endl;
	//cout << "ID of current_population_datas_: " << typeid(current_population_datas_).name() << endl;


	// ------------------------------------------------------------------------------------------------------------------
	*/
	for (ReporterPtrs::iterator reporter=config_.reporters.begin(); reporter!=config_.reporters.end(); ++reporter)
	{
		const bool is_final_generation = false;
		(*reporter)->update(current_generation_index_, *current_populations_, *current_population_datas_, is_final_generation);
	}

	if (os_popconfigs_)
	{
		os_popconfigs_ << "generation " << current_generation_index_ << endl;
		copy(popconfigs.begin(), popconfigs.end(), 
			 ostream_iterator<Population::Config>(os_popconfigs_,"\n"));
		os_popconfigs_ << endl;
	}

	++current_generation_index_;
}


void Simulator::simulate_all()
{
	if (config_.write_popconfig)
	{
		os_popconfigs_.open(bfs::path(config_.output_directory) / "forqs.popconfig.txt");
		if (!os_popconfigs_)
			throw runtime_error("[Simulator] Unable to open forqs.popconfig.txt");
	}

	const size_t generation_count = config_.population_config_generator->generation_count();
	for (size_t generation=0; generation <= generation_count; ++generation)
	{
		simulate_single_generation();
	}

	if (os_popconfigs_)
		os_popconfigs_.close();
}


void Simulator::update_final()
{
	for (ReporterPtrs::iterator reporter=config_.reporters.begin(); reporter!=config_.reporters.end(); ++reporter)
	{
		const bool is_final_generation = true;
		(*reporter)->update(current_generation_index_, *current_populations_, *current_population_datas_, is_final_generation);
	}
}


