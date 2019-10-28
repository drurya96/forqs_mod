//
// ReporterImplementation.hpp
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


#ifndef _REPORTERIMPLEMENTATION_HPP_
#define _REPORTERIMPLEMENTATION_HPP_


#include "Reporter.hpp"
#include "Simulator.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <ctime>


namespace bfs = boost::filesystem;


///
/// \defgroup Reporters Reporters
///
/// classes for reporting output from the simulation
///
/// Reporters are called every generation with data from the simulation.  Each Reporter
/// produces output files with filenames prefixed by the name of the Reporter.
///


//
// Reporter_Timer 
//

///
/// reports generation times
/// 
/// Parameters: none
///
/// Example: [example_stepping_stone.txt](../../examples/example_stepping_stone.txt)
///
/// \ingroup Reporters
///

class Reporter_Timer : public Reporter
{
    public:

    Reporter_Timer(const std::string& id);

    const std::vector<double>& times() {return times_;}
    double mean_generation_time() const;

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_Timer";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    
    private:

    std::clock_t begin_;
    std::vector<double> times_;
};



//
// Reporter_Population
//

///
/// reports full data for each population in each generation
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// update_step = \<int\> | 0 | optional, number of generations between updates (0 == final generation only)
///
/// Example: [example_neutral_admixture.txt](../../examples/example_neutral_admixture.txt)
///
/// \ingroup Reporters
///

class Reporter_Population : public Reporter
{
    public:

    Reporter_Population(const std::string& id);

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_Population";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    size_t update_step_;
};


//
// Reporter_AlleleFrequencies
//

///
/// reports allele frequencies for specified loci for each generation
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus = \<id_locus\> | none | required
/// locus_list = \<id_locus_list\> | none | required
/// quantitative_trait = \<id_quantitative_trait\> | none | required
///
/// Notes: 
///  - loci can be specified individually (Locus), as a LocusList, and/or by
///    QuantitativeTrait (all QTLs for the trait)
///  - assumes biallelic SNPs (variant values are 0/1)
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup Reporters
///

class Reporter_AlleleFrequencies : public Reporter
{
    // note: this handles binary SNPs only (0/1 values)

    public:

    Reporter_AlleleFrequencies(const std::string& id);

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    virtual Loci loci(size_t generation_index, 
                      bool is_final_generation) const {return loci_;}

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_AlleleFrequencies";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    Loci loci_;

    Loci loci_specified_;
    LocusListPtrs locus_lists_specified_;
    QuantitativeTraitPtrs qts_specified_;

    typedef shared_ptr<std::ostream> OstreamPtr;
    typedef std::map<Locus, OstreamPtr> OstreamMap;
    OstreamMap os_map_;

    std::string generate_filename(const Locus& locus) const;
    void open_streams();

    bool report_D_;
};


//
// Reporter_LD
//

///
/// reports linkage disequilibrium for a specified locus for each generation
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus_1 = \<id\> | none | required
/// locus_2 = \<id\> | none | required
///
/// Example: [example_ld.txt](../../examples/example_ld.txt)
///
/// \ingroup Reporters
///

class Reporter_LD : public Reporter
{
    // note: this handles binary SNPs only (0/1 values)

    public:

    Reporter_LD(const std::string& id, 
                Locus locus_1 = Locus("id_dummy_1"),
                Locus locus_2 = Locus("id_dummy_2"));

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    virtual Loci loci(size_t generation_index, 
                      bool is_final_generation) const;

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_LD";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    Locus locus_1_;
    Locus locus_2_;

    bfs::ofstream os_;

    void open_streams();
};


//
// Reporter_TraitValues
//

///
/// reports mean trait values for a specified quantitative trait,
/// for each population in each generation
///
/// parameter | default | notes
/// ----------|---------|-------------
/// quantitative_traits = \<id\> [...] | none | trait ids to report (required)
/// ignore_zero_values = \<id\> [...] | none | for trait ids in this list, mean will be computed on nonzero values only
/// write_full = 1 | 0 | write individual trait value data (per generation, per population)
/// filetag = \<string\> | tag added to filenames for multiple Reporter_TraitValues instances
///
/// Example: [example_stepping_stone.txt](../../examples/example_stepping_stone.txt)
///
/// \ingroup Reporters
///

class Reporter_TraitValues : public Reporter
{
    public:

    Reporter_TraitValues(const std::string& id);

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_TraitValues";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::vector<std::string> qtids_;
    std::set<std::string> qtids_ignore_zero_values_;
    bool write_full_;
    std::string filetag_;

    typedef shared_ptr<std::ostream> OstreamPtr;
    typedef std::vector<OstreamPtr> OstreamPtrs;
    OstreamPtrs os_means_;

    void open_streams();
};


//
// Reporter_HaplotypeDiversity
//

// in development
//
// reports haplotype diversity for each population in each generation
//
// Parameters:
//     - chromosome = \<chromosome_pair_index\> \<length\> \<step\>
//
// \ingroup Reporters
//

class Reporter_HaplotypeDiversity : public Reporter
{
    public:

    struct ChromosomeEntry
    {
        size_t index;
        size_t length;
        size_t step;

        ChromosomeEntry(size_t i = 0, size_t l = 0, size_t s = 0);
        ChromosomeEntry(const std::string& configuration);
        std::string configuration() const;
    };

    typedef std::map<size_t, ChromosomeEntry> ChromosomeEntries; // chromosome index -> ChromosomeEntry

    Reporter_HaplotypeDiversity(const std::string& id, 
                                const ChromosomeEntries& chromosome_entries = ChromosomeEntries());

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);
    
    // Configurable interface

    virtual std::string class_name() const {return "Reporter_HaplotypeDiversity";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    ChromosomeEntries chromosome_entries_;

    // map chromosome pair index -> vector<stream> (one stream for each population)
    typedef std::map< size_t, std::vector< shared_ptr<bfs::ofstream> > > ChromosomeStreams;
    ChromosomeStreams chromosome_streams_;

    void open_streams(size_t population_count);
};



//
// HaplotypeGrouping
//


class HaplotypeGrouping : public Configurable
{
    public:

    HaplotypeGrouping(const std::string& id) : Configurable(id) {}

    virtual size_t group_count() const = 0;
    virtual size_t group(unsigned int chromosome_id) const = 0;
    virtual ~HaplotypeGrouping() {}

    // Configurable interface

    virtual std::string class_name() const;
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
};


typedef shared_ptr<HaplotypeGrouping> HaplotypeGroupingPtr;


//
// HaplotypeGrouping_IDRange
//

///
/// basic grouping by haplotype id, for use with Reporter_HaplotypeFrequencies
///
/// parameter | default | notes
/// ----------|---------|-------------
/// start:count = \<int_id_start\> \<int_id_count\> | none | multiple allowed
///
/// Example: [example_2_locus_selection.txt](../../examples/example_2_locus_selection.txt)
///
/// \ingroup Reporters
///


class HaplotypeGrouping_IDRange : public HaplotypeGrouping
{
    public:

    HaplotypeGrouping_IDRange(const std::string& id);

    virtual size_t group_count() const;
    virtual size_t group(unsigned int chromosome_id) const;

    // Configurable interface

    virtual std::string class_name() const {return "HaplotypeGrouping_IDRange";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);

    private:

    std::vector<unsigned int> id_starts_;
    std::vector<unsigned int> id_counts_;
};


//
// HaplotypeGrouping_Uniform
//

///
/// haplotype groups defined by the number of ids per group
///
/// parameter | default | notes
/// ----------|---------|-------------
/// ids_per_group = \<int\> | none | population size must be a multiple of this value
///
/// Note: single population only
///
/// Example: [tutorial_6_haplotype_frequencies.txt](../../examples/tutorial_6_haplotype_frequencies.txt)
///
/// \ingroup Reporters
///


class HaplotypeGrouping_Uniform : public HaplotypeGrouping
{
    public:

    HaplotypeGrouping_Uniform(const std::string& id);

    virtual size_t group_count() const;
    virtual size_t group(unsigned int chromosome_id) const;

    // Configurable interface

    virtual std::string class_name() const {return "HaplotypeGrouping_Uniform";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);

    private:
    unsigned int ids_per_group_; 
    unsigned int id_offset_;
    unsigned int id_count_;
};


//
// Reporter_HaplotypeFrequencies
//

///
/// reports local haplotype frequencies for each population
///
/// parameter | default | notes
/// ----------|---------|-------------
/// haplotype_grouping = \<id\> | none | required
/// chromosome_step = \<int\> | none | required
/// update_step = \<int\> | 0 | optional, number of generations between updates (0 == final generation only)
///
/// Example: [example_2_locus_selection.txt](../../examples/example_2_locus_selection.txt)
///
/// \ingroup Reporters
///


class Reporter_HaplotypeFrequencies : public Reporter
{
    public:

    Reporter_HaplotypeFrequencies(const std::string& id, 
                                  HaplotypeGroupingPtr haplotype_grouping = HaplotypeGroupingPtr(),
                                  size_t chromosome_step = 0,
                                  size_t update_step = 0);

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);
   
    // Configurable interface

    virtual std::string class_name() const {return "Reporter_HaplotypeFrequencies";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    HaplotypeGroupingPtr haplotype_grouping_;
    size_t chromosome_step_;
    size_t update_step_;
    std::vector<unsigned int> chromosome_lengths_;

    void write_haplotype_frequencies(std::ostream& os, const Population& population,
        size_t chromosome_pair_index, unsigned int chromosome_length) const;
};


//
// Reporter_Regions
//

///
/// reports haplotypes within a region
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus:length = \<id_locus_start\> \<int_region_length\> | none | multiple allowed
///
/// Example: [example_neutral_mutation_region.txt](../../examples/example_neutral_mutation_region.txt)
///
/// \ingroup Reporters
///


class Reporter_Regions : public Reporter
{
    public:

    struct Region
    {
        Locus locus;
        size_t length;

        Region(const Locus& _locus = Locus("dummy_id"),
               size_t _length = 0)
        :   locus(_locus), length(_length)
        {}

        std::string configuration() const;
    };

    typedef std::vector<Region> Regions;

    Reporter_Regions(const std::string& id, 
                     const Regions& regions = Regions());

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);
   
    // Configurable interface

    virtual std::string class_name() const {return "Reporter_Regions";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    private:

    Regions regions_;
    unsigned int ms_mapping_begin_;
    unsigned int ms_mapping_end_;
};


//
// Reporter_DeterministicTrajectories
//

///
/// convenience reporter outputs deterministic allele frequencies for specified genotypic fitness values
///
/// parameter | default | notes
/// ----------|---------|-------------
/// initial_allele_frequency = \<float\> | none | required
/// w0 = \<float\> | none | required
/// w1 = \<float\> | none | required
/// w2 = \<float\> | none | required
///
/// Alternative parametrization:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// initial_allele_frequency = \<float\> | none | required
/// additive_selection_coefficient = \<float\> | none | required; (w0,w1,w2) == (1,1+s,1+2s)
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup Reporters
///


class Reporter_DeterministicTrajectories : public Reporter
{
    public:

    Reporter_DeterministicTrajectories(const std::string& id,
                                       double initial_allele_frequency = 0,
                                       std::vector<double> w = std::vector<double>(3))
    :   Configurable(id), initial_allele_frequency_(initial_allele_frequency), w_(w), generation_count_(0)
    {}

    virtual void update(size_t generation_number,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_DeterministicTrajectories";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    
    private:
    
    double initial_allele_frequency_;
    std::vector<double> w_;
    double generation_count_;
};


//
// Reporter_Variants 
//

///
/// reports variant values for all individuals
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// loci = \<id\> [...] | none | list of ids (Locus or LocusList)
///
/// Example: [turner_minimal.txt](../../examples/turner_minimal.txt)
///
/// \ingroup Reporters
///

class Reporter_Variants : public Reporter
{
    public:

    Reporter_Variants(const std::string& id);

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    virtual Loci loci(size_t generation_index, bool is_final_generation) const;

    // Configurable interface

    virtual std::string class_name() const {return "Reporter_Variants";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;
    
    private:

    std::vector<std::string> locus_ids_;
    size_t update_step_;
    std::vector<Locus> loci_;
    ConfigurablePtrs children_;
    QuantitativeTraitPtrs qts_specified_;
};


#endif // _REPORTERIMPLEMENTATION_HPP_

