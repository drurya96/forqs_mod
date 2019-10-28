//
// VariantIndicatorImplementation.hpp
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


#ifndef _VARIANTINDICATORIMPLEMENTATION_HPP_
#define _VARIANTINDICATORIMPLEMENTATION_HPP_


#include "VariantIndicator.hpp"
#include "MSFormat.hpp"
#include "Random.hpp"
#include "Reporter.hpp"


///
/// \defgroup VariantIndicators VariantIndicators
///
/// classes for specifying haplotype SNP values
///


//
// VariantIndicator_Trivial
//

///
/// trivial implementation: returns 0
/// 
/// Parameters: none
///
/// \ingroup VariantIndicators
///

class VariantIndicator_Trivial : public VariantIndicator
{
    public:

    VariantIndicator_Trivial(const std::string& id) : Configurable(id) {} 

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const {return 0;}

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_Trivial";}
    virtual Parameters parameters() const {return Parameters();}
    virtual void configure(const Parameters& parameters, const Registry& registry) {}
};


//
// VariantIndicator_Composite
//

///
/// encapsulates a list of VariantIndicators; returns first non-zero value for
/// a given (haplotype id,locus) if one exists
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// variant_indicators = \<id\> | none | children are called in order specified
///
/// Example: TODO
///
/// \ingroup VariantIndicators
///

class VariantIndicator_Composite : public VariantIndicator
{
    public:

    VariantIndicator_Composite(const std::string& id) : Configurable(id) {} 

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_Composite";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    private:

    VariantIndicatorPtrs variant_indicators_;
};


//
// VariantIndicator_IDRange
//

///
/// returns variant values based on specified haplotype id ranges and loci
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// locus:start:count:step:value = \<id_locus\> \<int_id_start\> \<int_id_count\> \<int_id_step\> \<int_value\> | none | multiple allowed
///
/// Example: [example_2_locus_selection.txt](../../examples/example_2_locus_selection.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_IDRange : public VariantIndicator
{
    public:

    VariantIndicator_IDRange(const std::string& id) : Configurable(id) {} 

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_IDRange";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    protected:

    struct Entry
    {
        unsigned int id_start;
        unsigned int id_count;
        unsigned int id_step;
        unsigned int value;

        Entry(unsigned int _id_start, unsigned int _id_count, unsigned int _id_step, unsigned int _value)
        :   id_start(_id_start), id_count(_id_count), id_step(_id_step), value(_value)
        {}
    };

    typedef std::multimap<Locus,Entry> EntryMap;
    EntryMap entries_;
};


//
// VariantIndicator_IDSet
//

///
/// flexible implementation:  user specifies per-locus value and id list
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// locus:value:ids = \<id_locus\> \<int_value\> \<int_id\> [\<int_id\> ...] | none | multiple allowed
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_IDSet : public VariantIndicator
{
    public:

    VariantIndicator_IDSet(const std::string& id) : Configurable(id) {} 
    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;
    virtual void write_file(const std::string& filename) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_IDSet";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    struct Entry
    {
        unsigned int value;
        std::set<unsigned int> ids;

        Entry(unsigned int _value = -1u, std::set<unsigned int> _ids = std::set<unsigned int>())
        :   value(_value), ids(_ids)
        {}
    };

    typedef std::map<Locus,Entry> EntryMap;

    const EntryMap& entries() const {return entries_;} // for testing

    protected:

    EntryMap entries_;
};


//
// VariantIndicator_Random
//

///
/// Assigns SNP values randomly to individuals in the initial populations,
/// independently for each locus with specified allele frequencies,
/// so that each SNP will be in Hardy-Weinberg equilibrium (in expectation).
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus_list:population:frequencies = \<id_locus\> \<int_population\> \<float_frequency\> [...] | none | multiple allowed; number of frequencies must match number of loci
/// locus_list:population:frequency_distribution = \<id_locus\> \<int_population\> \<id_distribution\> [...] | none | multiple allowed
///
/// Note: * may be used to indicate all populations
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_Random : public VariantIndicator_IDSet
{
    public:

    VariantIndicator_Random(const std::string& id) : Configurable(id), VariantIndicator_IDSet(id) {}

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_Random";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& simconfig);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    private:

    struct LocusListInfo
    {
        size_t population_index;
        LocusListPtr locus_list;
        std::vector<double> frequencies;
        Random::DistributionPtr distribution;

        LocusListInfo(size_t _population_index, const LocusListPtr& _locus_list, 
                      std::vector<double> _frequencies)
        :   population_index(_population_index), locus_list(_locus_list), frequencies(_frequencies) 
        {}

        LocusListInfo(size_t _population_index, const LocusListPtr& _locus_list, 
                      Random::DistributionPtr _distribution)
        :   population_index(_population_index), locus_list(_locus_list), distribution(_distribution) 
        {}
    };

    typedef std::vector<LocusListInfo> LocusListInfos;
    LocusListInfos locus_list_infos_;
};


//
// VariantIndicator_File
//

///
/// returns variant values based on an ms-format file
/// 
/// parameter | default | notes
/// ----------|---------|-------------
/// msfile = \<string_filename\> | none | required
/// loci = \<id\> ... | none | multiple Locus or LocusList ids allowed
///
/// Note:  The 'positions:' values defined in the ms-format file are ignored; the 
///        loci specified as parameters (loci = \<id\> ...) are used instead.
///
/// Example: [example_1_locus_selection_msfile.txt](../../examples/example_1_locus_selection_msfile.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_File : public VariantIndicator
{
    public:

    VariantIndicator_File(const std::string& id) : Configurable(id) {} 

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_File";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    private:

    std::string msfile_;
    MSFormatPtr ms_;
    std::vector<std::string> locus_ids_;

    std::vector<Locus> loci_;
    std::map<Locus, size_t> locus_index_map_;
    ConfigurablePtrs children_;
};


//
// VariantIndicator_SingleLocusHardyWeinberg
//

///
/// Assigns SNP value 1 to both haplotypes of the first Np^2 individuals, and to
/// one haplotype of the next N*2pq individuals. 
///
/// N == population_size, p == 1-q == allele_frequency
///
/// Assumes chromosome ids are assigned consecutively (2 per individual): 0, ..., 2N-1
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus = \<id\> | none | required
/// allele_frequency = \<float\> | none | required
///
/// Implementation note: requires call to initialize() to obtain population id ranges from
///                      PopulationConfigGenerator
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_SingleLocusHardyWeinberg : public VariantIndicator_IDRange
{
    public:

    VariantIndicator_SingleLocusHardyWeinberg(const std::string& id, 
                                              Locus locus = Locus("id_dummy"),
                                              double allele_frequency = 0);
    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_SingleLocusHardyWeinberg";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& simconfig);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    private:

    Locus locus_;
    double allele_frequency_;
};


//
// VariantIndicator_TwoLocusLD
//

///
/// Assigns SNP values according to allele frequencies at two loci and specified
/// linkage disequilibrium (D).  Note:  all individuals are homozygotes, i.e. their
/// haplotypes are identical.
///
/// Assumes chromosome ids are assigned consecutively (2 per individual): 0, ..., 2N-1
///
/// parameter | default | notes
/// ----------|---------|-------------
/// population_size = \<int\> | none | required
/// locus_1 = \<id\> | none | required
/// allele_frequency_1 = \<float\> | none | required
/// locus_2 = \<id\> | none | required
/// allele_frequency_2 = \<float\> | none | required
/// id_offset_step = \<int\>) | none | required
///
/// Example: [example_ld.txt](../../examples/example_ld.txt)
///
/// \ingroup VariantIndicators
///

class VariantIndicator_TwoLocusLD : public VariantIndicator
{
    public:

    VariantIndicator_TwoLocusLD(const std::string& id, 
                            size_t population_size = 0,
                            Locus locus_1 = Locus("id_dummy_1"),
                            double allele_frequency_1 = 0,
                            Locus locus_2 = Locus("id_dummy_2"),
                            double allele_frequency_2 = 0,
                            double D = 0,
                            unsigned int id_offset_step = 0);

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_TwoLocusLD";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    private:

    size_t population_size_;
    Locus locus_1_;
    double allele_frequency_1_;
    Locus locus_2_;
    double allele_frequency_2_;
    double D_;
    unsigned int id_offset_step_;

    size_t max_00_;
    size_t max_01_;
    size_t max_10_;
    size_t max_11_;

    void initialize_internal();
};


//
// VariantIndicator_Mutable
//

//
// wrapper implementation for handling mutation
// 

class VariantIndicator_Mutable : public VariantIndicator, public Reporter
{
    public:

    VariantIndicator_Mutable(const std::string& id,
                             unsigned int unused_id_start = 0,
                             VariantIndicatorPtr vi = VariantIndicatorPtr(),
                             const std::string& output_directory = "");

    virtual unsigned int operator()(unsigned int chunk_id, const Locus& locus) const;

    virtual unsigned int mutate(unsigned int old_chunk_id, const Locus& locus, unsigned int value);

    void report(std::ostream& os) const;

    // Configurable interface

    virtual std::string class_name() const {return "VariantIndicator_Mutable";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& simconfig);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& written) const;

    // Reporter interface

    virtual void update(size_t generation_index,
                        const PopulationPtrs& populations,
                        const PopulationDataPtrs& population_datas,
                        bool is_final_generation);

    virtual Loci loci(size_t generation_index, 
                      bool is_final_generation) const;

    private:

    unsigned int unused_id_start_;
    unsigned int unused_id_current_;
    VariantIndicatorPtr vi_;
    bfs::path outdir_;
    bfs::ofstream os_debug_;

    typedef std::map<unsigned int, unsigned int> IDAncestry;
    IDAncestry id_ancestry_;

    typedef std::map<unsigned int, unsigned int> IDValueMap;
    typedef std::map<Locus, IDValueMap> IDValueMaps;
    IDValueMaps id_value_maps_;
};


#endif // _VARIANTINDICATORIMPLEMENTATION_HPP_


