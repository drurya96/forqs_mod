//
// ReporterImplementation.cpp
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


#include "ReporterImplementation.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <stdexcept>
#include <numeric>
#include <sstream>


using namespace std;
namespace bfs = boost::filesystem;


//
// Reporter_Timer
//


Reporter_Timer::Reporter_Timer(const string& id)
:   Configurable(id), begin_(clock())
{}


void Reporter_Timer::update(size_t generation_index,
                            const PopulationPtrs& populations,
                            const PopulationDataPtrs& population_datas,
                            bool is_final_generation)
{
    times_.push_back((clock() - begin_)/double(CLOCKS_PER_SEC));

    if (is_final_generation)
    {
        if (!output_directory_.empty())
        {
            bfs::ofstream os_(output_directory_ / "timer.txt");
            if (!os_) throw runtime_error("[Reporter_Timer] Unable to open timer.txt");
            copy(times_.begin(), times_.end(), ostream_iterator<double>(os_, "\n"));
            os_ << "mean_generation_time: " << mean_generation_time() << endl;
            os_.close();
        }
    }
}


double Reporter_Timer::mean_generation_time() const
{
    // report the mean time between successive update() calls;
    // ignore first (times_[0]) and final update() call (times_[size-1])

    if (times_.size() < 2)
        throw runtime_error("[Reporter_Timer::mean_generation_time()] No times.");

    vector<double> diffs;
    adjacent_difference(times_.begin(), times_.end(), back_inserter(diffs));
    return accumulate(diffs.begin()+1, diffs.end()-1, 0.0) / (diffs.size()-2);
}


Parameters Reporter_Timer::parameters() const
{
    return Parameters();
}


void Reporter_Timer::configure(const Parameters& parameters, const Registry& registry)
{
}


//
// Reporter_Population
//


Reporter_Population::Reporter_Population(const string& id)
:   Configurable(id), update_step_(0)
{}


void Reporter_Population::update(size_t generation_index,
                                 const PopulationPtrs& populations,
                                 const PopulationDataPtrs& population_datas,
                                 bool is_final_generation)
{
    if (update_step_ && (generation_index % update_step_ == 0) || is_final_generation)
    {
        if (populations.size() != population_datas.size())
            throw runtime_error("[Reporter_Population] Population data size mismatch.");

        const size_t population_count = populations.size();
        const char* filestem = "population";

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;

            if (is_final_generation)
                filename << filestem << "_final_pop" << population_index + 1 << ".txt"; 
            else
                filename << filestem << "_gen" << generation_index << "_pop" << population_index + 1 << ".txt"; 

            bfs::ofstream os(output_directory_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_Population] Unable to open " + filename.str()).c_str());

            os << *populations[population_index];
            os.close();
        }
    }
}


Parameters Reporter_Population::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("update_step", update_step_);
    return parameters;
}


void Reporter_Population::configure(const Parameters& parameters, const Registry& registry)
{
    update_step_ = parameters.value<size_t>("update_step", 0);
}


//
// Reporter_AlleleFrequencies
//


Reporter_AlleleFrequencies::Reporter_AlleleFrequencies(const string& id)
:   Configurable(id), report_D_(false)
{}


void Reporter_AlleleFrequencies::update(size_t generation_index,
                                        const PopulationPtrs& populations,
                                        const PopulationDataPtrs& population_datas,
                                        bool is_final_generation)
{
    if (is_final_generation)
    {
        bfs::ofstream os(output_directory_ / (object_id() + "_summary.txt"));
        if (!os)
            throw runtime_error("[Reporter_AlleleFrequencies] Unable to open allele_frequencies_final_summary.txt");

        os << "chromosome position ";
        for (size_t i=0; i<population_datas.size(); ++i)
            os << "pop" << i << " ";
        if (report_D_ && population_datas.size() == 4) os << "D";
        os << endl;

        for (Loci::const_iterator locus=loci_.begin(); locus!=loci_.end(); ++locus)
        {
            os << locus->chromosome_pair_index + 1 << " " << locus->position << " ";

            vector<double> allele_freqs;

            for (PopulationDataPtrs::const_iterator data=population_datas.begin(); data!=population_datas.end(); ++data)
            {
                GenotypeDataPtr genotypes = (*data)->genotypes->get(*locus);
                allele_freqs.push_back(genotypes->allele_frequency());
            }

            copy(allele_freqs.begin(), allele_freqs.end(), ostream_iterator<double>(os," "));

            if (report_D_ && allele_freqs.size() == 4)
            {
                double D = allele_freqs[2] - allele_freqs[0] + allele_freqs[3] - allele_freqs[1];
                os << D;
            }

            os << endl;
        }

        return;
    }

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_AlleleFrequencies] Population data size mismatch.");

    if (os_map_.empty()) open_streams();

    for (Loci::const_iterator locus=loci_.begin(); locus!=loci_.end(); ++locus)
    {
        ostream& os = *os_map_[*locus];
        for (PopulationDataPtrs::const_iterator data=population_datas.begin(); data!=population_datas.end(); ++data)
        {
            GenotypeDataPtr genotypes = (*data)->genotypes->get(*locus);
            os << genotypes->allele_frequency() << " "; // note: this is for binary allele case only
        }
        os << endl;
    }
}


Parameters Reporter_AlleleFrequencies::parameters() const
{
    Parameters parameters;
    for (Loci::const_iterator it=loci_specified_.begin(); it!=loci_specified_.end(); ++it)
        parameters.insert_name_value("locus", it->object_id());
    for (LocusListPtrs::const_iterator it=locus_lists_specified_.begin(); it!=locus_lists_specified_.end(); ++it)
        parameters.insert_name_value("locus_list", (*it)->object_id());
    for (QuantitativeTraitPtrs::const_iterator it=qts_specified_.begin(); it!=qts_specified_.end(); ++it)
        parameters.insert_name_value("quantitative_trait", (*it)->object_id());
    return parameters;
}


void Reporter_AlleleFrequencies::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> locus_ids = parameters.values<string>("locus");

    for (vector<string>::const_iterator id=locus_ids.begin(); id!=locus_ids.end(); ++id)
    {
        const LocusPtr& locus = registry.get<Locus>(*id);
        loci_specified_.insert(*locus);
    }

    vector<string> locus_list_ids = parameters.values<string>("locus_list");

    for (vector<string>::const_iterator llid=locus_list_ids.begin(); llid!=locus_list_ids.end(); ++llid)
    {
        const LocusListPtr& locus_list = registry.get<LocusList>(*llid);
        locus_lists_specified_.push_back(locus_list);
    }

    vector<string> qt_ids = parameters.values<string>("quantitative_trait");

    for (vector<string>::const_iterator qtid=qt_ids.begin(); qtid!=qt_ids.end(); ++qtid)
    {
        const QuantitativeTraitPtr& qt = registry.get<QuantitativeTrait>(*qtid);
        qts_specified_.push_back(qt);
    }

    report_D_ = parameters.value<bool>("report_D", false);
}


void Reporter_AlleleFrequencies::initialize(const SimulatorConfig& config)
{
    Reporter::initialize(config);

    loci_.clear();
    for (Loci::const_iterator it=loci_specified_.begin(); it!=loci_specified_.end(); ++it)
        loci_.insert(*it);
    for (LocusListPtrs::const_iterator it=locus_lists_specified_.begin(); it!=locus_lists_specified_.end(); ++it)
        copy((*it)->begin(), (*it)->end(), inserter(loci_, loci_.end()));
    for (QuantitativeTraitPtrs::const_iterator it=qts_specified_.begin(); it!=qts_specified_.end(); ++it)
        copy((*it)->loci().begin(), (*it)->loci().end(), inserter(loci_, loci_.end()));
}


void Reporter_AlleleFrequencies::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (Loci::const_iterator it=loci_specified_.begin(); it!=loci_specified_.end(); ++it)
        it->write_configuration(os, ids_written);
    for (LocusListPtrs::const_iterator it=locus_lists_specified_.begin(); it!=locus_lists_specified_.end(); ++it)
        (*it)->write_configuration(os, ids_written);
    for (QuantitativeTraitPtrs::const_iterator it=qts_specified_.begin(); it!=qts_specified_.end(); ++it)
        (*it)->write_configuration(os, ids_written);
}


string Reporter_AlleleFrequencies::generate_filename(const Locus& locus) const
{
    ostringstream filename;
    filename << "allele_frequencies"
             << "_chr" << locus.chromosome_pair_index + 1
             << "_pos" << locus.position 
             << ".txt";
    return filename.str();
}


void Reporter_AlleleFrequencies::open_streams()
{
    for (Loci::const_iterator locus=loci_.begin(); locus!=loci_.end(); ++locus)
    {
        bfs::path filename = output_directory_ / generate_filename(*locus);
        OstreamPtr os(new bfs::ofstream(filename));
        if (!os)
            throw runtime_error(("[Reporter_AlleleFrequencies] Unable to open file " + filename.string()).c_str());
        os_map_[*locus] = os;
    }
}


//
// Reporter_LD
//


Reporter_LD::Reporter_LD(const string& id,
                         Locus locus_1,
                         Locus locus_2)
:   Configurable(id), locus_1_(locus_1), locus_2_(locus_2)
{}


void Reporter_LD::update(size_t generation_index,
                         const PopulationPtrs& populations,
                         const PopulationDataPtrs& population_datas,
                         bool is_final_generation)
{
    if (is_final_generation) return;

    if (!os_.is_open()) open_streams();

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_LD] Population data size mismatch.");

    for (PopulationDataPtrs::const_iterator data=population_datas.begin(); data!=population_datas.end(); ++data)
    {
        const GenotypeData& genotypes1 = *(*data)->genotypes->get(locus_1_);
        const GenotypeData& genotypes2 = *(*data)->genotypes->get(locus_2_);

        if (genotypes1.size() != genotypes2.size())
            throw runtime_error("[Reporter_LD] Genotype vector size mismatch.");

        double counts[2][2] = {{0, 0}, {0, 0}};

        for (GenotypeData::const_iterator gt1=genotypes1.begin(), gt2=genotypes2.begin();
             gt1!=genotypes1.end(); ++gt1, ++gt2)
        {
            ++counts[(size_t)genotype_first(*gt1)][(size_t)genotype_first(*gt2)];
            ++counts[(size_t)genotype_second(*gt1)][(size_t)genotype_second(*gt2)];
        }

        double count_total = counts[0][0] + counts[0][1] + counts[1][0] + counts[1][1];
        double D = (counts[0][0]*counts[1][1] - counts[0][1]*counts[1][0])/count_total/count_total;

        /*
        os_ << "(" << counts[0][0]/count_total << ","
                   << counts[0][1]/count_total << ","
                   << counts[1][0]/count_total << ","
                   << counts[1][1]/count_total << ") ";
        */

        os_ << D << " ";
    }
    os_ << endl;
}


Loci Reporter_LD::loci(size_t generation_index, bool is_final_generation) const
{
    Loci loci;
    loci.insert(locus_1_);
    loci.insert(locus_2_);
    return loci;
}


Parameters Reporter_LD::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("locus_1", locus_1_.object_id());
    parameters.insert_name_value("locus_2", locus_2_.object_id());
    return parameters;
}


void Reporter_LD::configure(const Parameters& parameters, const Registry& registry)
{
    locus_1_ = *registry.get<Locus>(parameters.value<string>("locus_1"));
    locus_2_ = *registry.get<Locus>(parameters.value<string>("locus_2"));
}


void Reporter_LD::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    locus_1_.write_configuration(os, ids_written);
    locus_2_.write_configuration(os, ids_written);
}


void Reporter_LD::open_streams()
{
    ostringstream filename;
    filename << "ld"
             << "_chr" << locus_1_.chromosome_pair_index + 1
             << "_pos" << locus_1_.position 
             << "_chr" << locus_2_.chromosome_pair_index + 1
             << "_pos" << locus_2_.position 
             << ".txt";

    os_.open(output_directory_ / filename.str());
    if (!os_)
        throw runtime_error(("[Reporter_LD] Unable to open file " + filename.str()).c_str());
}


//
// Reporter_TraitValues
//


Reporter_TraitValues::Reporter_TraitValues(const string& id)
:   Configurable(id), write_full_(false)
{}


void Reporter_TraitValues::update(size_t generation_index,
                                  const PopulationPtrs& populations,
                                  const PopulationDataPtrs& population_datas,
                                  bool is_final_generation)
{
    if (is_final_generation) return;

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_TraitValues] Population data size mismatch.");

    if (os_means_.empty()) open_streams();

    // one file per trait, containing mean trait values for all generations & populations

    OstreamPtrs::iterator os = os_means_.begin();
    for (vector<string>::const_iterator id=qtids_.begin(); id!=qtids_.end(); ++id, ++os)
    {
        for (PopulationDataPtrs::const_iterator data=population_datas.begin(); data!=population_datas.end(); ++data)
        {
            const DataVector& trait_values = *(*data)->trait_values->get(*id);

            // update os_mean_

            if (!qtids_ignore_zero_values_.count(*id))
            {
                **os << trait_values.mean() << " ";
            }
            else
            {
                DataVector nonzero_trait_values;
                remove_copy_if(trait_values.begin(), trait_values.end(), 
                    back_inserter(nonzero_trait_values), bind2nd(equal_to<double>(), 0.));
                **os << nonzero_trait_values.mean() << " ";
            }

        }
        **os << endl;
    }

    if (write_full_)
    {
        // one file per generation per population, containing individual values for all traits

        for (PopulationDataPtrs::const_iterator data=population_datas.begin(); data!=population_datas.end(); ++data)
        {
            const size_t population_number = data - population_datas.begin() + 1; // 1-based
            const size_t population_size = (*data)->population_size;

            ostringstream filename;
            filename << "trait_values_";
            if (!filetag_.empty()) filename << filetag_ << "_";
            filename << "full_gen" << generation_index
                     << "_pop" << population_number << ".txt";

            bfs::ofstream os(output_directory_ / filename.str());
            if (!os)
                throw runtime_error("[Reporter_TraitValues] Unable to open file " + filename.str());

            // iterate through all traits' DataVectors in parallel (output columns)

            vector<DataVector::const_iterator> its;

            for (vector<string>::const_iterator id=qtids_.begin(); id!=qtids_.end(); ++id)
            {
                const DataVector& trait_values = *(*data)->trait_values->get(*id);
                if (trait_values.size() != population_size)
                    throw runtime_error("[Reporter_TraitValues] Data size mismatch.");
                its.push_back(trait_values.begin());
            }

            copy(qtids_.begin(), qtids_.end(), ostream_iterator<string>(os, " "));
            os << endl;
            
            for (size_t i=0; i<population_size; ++i)
            {
                for (vector<DataVector::const_iterator>::iterator it=its.begin(); 
                     it!=its.end(); ++it)
                {
                    os << **it << "\t";
                    ++(*it);
                }
                os << endl;
            }
        }
    }
}


Parameters Reporter_TraitValues::parameters() const
{
    Parameters parameters;

    parameters.insert_name_value_vector("quantitative_traits", qtids_);

    if (!qtids_ignore_zero_values_.empty()) 
    {
        vector<string> qtids_ignore;
        copy(qtids_ignore_zero_values_.begin(), qtids_ignore_zero_values_.end(),
             back_inserter(qtids_ignore));
        parameters.insert_name_value_vector("ignore_zero_values", qtids_ignore);
    }

    parameters.insert_name_value("write_full", write_full_);

    if (!filetag_.empty())
        parameters.insert_name_value("filetag", filetag_);

    return parameters;
}


void Reporter_TraitValues::configure(const Parameters& parameters, const Registry& registry)
{
    qtids_ = parameters.value_vector<string>("quantitative_traits");

    vector<string> qtids_ignore = parameters.value_vector<string>("ignore_zero_values", vector<string>());
    copy(qtids_ignore.begin(), qtids_ignore.end(), 
            inserter(qtids_ignore_zero_values_, qtids_ignore_zero_values_.begin()));

    write_full_ = parameters.value<int>("write_full", false);
    filetag_ = parameters.value<string>("filetag", "");
}


void Reporter_TraitValues::open_streams()
{
    for (vector<string>::const_iterator id=qtids_.begin(); id!=qtids_.end(); ++id)
    {
        ostringstream filename;
        filename << "trait_values_";
        if (!filetag_.empty()) filename << filetag_ << "_";
        filename << "mean_" << *id << ".txt";

        OstreamPtr os(new bfs::ofstream(output_directory_ / filename.str()));
        if (!*os)
            throw runtime_error("[Reporter_TraitValues] Unable to open file " + filename.str());

        os_means_.push_back(os);
    }
}


//
// Reporter_HaplotypeDiversity
//


Reporter_HaplotypeDiversity::ChromosomeEntry::ChromosomeEntry(size_t i, size_t l, size_t s)
:   index(i), length(l), step(s)
{}


Reporter_HaplotypeDiversity::ChromosomeEntry::ChromosomeEntry(const string& configuration)
{
    istringstream iss(configuration);
    iss >> index >> length >> step;
}


string Reporter_HaplotypeDiversity::ChromosomeEntry::configuration() const
{
    ostringstream oss;
    oss << index << " " << length << " " << step;
    return oss.str();
}


Reporter_HaplotypeDiversity::Reporter_HaplotypeDiversity(const string& id, 
                                                         const ChromosomeEntries& chromosome_entries)
:   Configurable(id), 
    chromosome_entries_(chromosome_entries)
{}


void Reporter_HaplotypeDiversity::open_streams(size_t population_count)
{
    if (output_directory_.string().empty())
        throw runtime_error("[Reporter_HaplotypeDiversity] No output directory specified.");

    // open a stream for each chromosome entry, for each population

    for (ChromosomeEntries::const_iterator it=chromosome_entries_.begin(); it!=chromosome_entries_.end(); ++it)
    {
        size_t chromosome_pair_index = it->first;

        chromosome_streams_[chromosome_pair_index] = vector< shared_ptr<bfs::ofstream> >();

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;
            filename << "haplotype_diversity_chr" << chromosome_pair_index + 1
                     << "_pop" << population_index + 1 << ".txt";

            shared_ptr<bfs::ofstream> os(new bfs::ofstream(output_directory_ / filename.str()));
            if (!*os)
                throw runtime_error(("[Reporter_HaplotypeDiversity] Unable to open file " + filename.str()).c_str());

            chromosome_streams_[chromosome_pair_index].push_back(os);

            const ChromosomeEntry& entry = it->second;

            *os << "# ";
            for (size_t position=0; position<entry.length; position+=entry.step)
                *os << position << " ";
            *os << endl;
        }
    }
}


void Reporter_HaplotypeDiversity::update(size_t generation_index,
                                         const PopulationPtrs& populations,
                                         const PopulationDataPtrs& population_datas,
                                         bool is_final_generation)
{
    if (is_final_generation) return;

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_HaplotypeDiversity] Population data size mismatch.");

    const size_t population_count = populations.size();

    if (chromosome_streams_.empty())
        open_streams(population_count);

    for (ChromosomeEntries::const_iterator it=chromosome_entries_.begin(); it!=chromosome_entries_.end(); ++it)
    {
        size_t chromosome_pair_index = it->first;
        const ChromosomeEntry& entry = it->second;

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            const Population& population = *populations[population_index];

            bfs::ofstream& os = *chromosome_streams_[chromosome_pair_index][population_index];

            //os << generation_index << " " << chromosome_pair_index + 1 << " " <<  population_index + 1  << endl; // TODO: remove

            for (size_t position=0; position<entry.length; position+=entry.step)
            {
                map<unsigned int, double> counts; // haplotype id -> count 
                double count_total = 0;

                ChromosomePairRangeIterator it = population.begin();
                ChromosomePairRangeIterator end = population.end();
                for (; it!=end; ++it)
                {
                    if (chromosome_pair_index >= it->size())
                        throw runtime_error("[Reporter_HaplotypeDiversity] Invalid chromosome pair index.");

                    const ChromosomePair& p = *(it->begin() + chromosome_pair_index);

                    unsigned int id1 = p.first.find_haplotype_chunk(position)->id;
                    unsigned int id2 = p.second.find_haplotype_chunk(position)->id;

                    ++counts[id1];
                    ++counts[id2];
                    count_total += 2;
                }

                os << counts.size() << " "; // for now, just report the number of different haplotypes
            }

            os << endl;
        }
    }
}


Parameters Reporter_HaplotypeDiversity::parameters() const
{
    Parameters parameters;
    for (ChromosomeEntries::const_iterator it=chromosome_entries_.begin(); it!=chromosome_entries_.end(); ++it)
        parameters.insert_name_value("chromosome", it->second.configuration());
    return parameters;
}


void Reporter_HaplotypeDiversity::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("chromosome");

    for (vector<string>::const_iterator configuration=configurations.begin(); configuration!=configurations.end(); ++configuration)
    {
        ChromosomeEntry entry(*configuration);
        chromosome_entries_[entry.index] = entry;
    }
}


//
// HaplotypeGrouping
//


string HaplotypeGrouping::class_name() const
{
    cerr << "[HaplotypeGrouping] Warning: virtual class_name() has not been defined in derived class.\n";
    return "HaplotypeGrouping";
}


Parameters HaplotypeGrouping::parameters() const
{
    cerr << "[HaplotypeGrouping] Warning: virtual parameters() has not been defined in derived class.\n";
    return Parameters();
}


void HaplotypeGrouping::configure(const Parameters& parameters, const Registry& registry)
{
    cerr << "[HaplotypeGrouping] Warning: virtual configure() has not been defined in derived class.\n";
}


//
// HaplotypeGrouping_IDRange
//


HaplotypeGrouping_IDRange::HaplotypeGrouping_IDRange(const string& id)
:   HaplotypeGrouping(id)
{}


size_t HaplotypeGrouping_IDRange::group_count() const
{
    return id_starts_.size();
}


size_t HaplotypeGrouping_IDRange::group(unsigned int chromosome_id) const
{
    if (id_starts_.empty())
        throw runtime_error("[HaplotypeGrouping_IDRange] No groups defined.");

    vector<unsigned int>::const_iterator it = upper_bound(id_starts_.begin(), id_starts_.end(), chromosome_id);

    if (it == id_starts_.begin())
    {
        ostringstream oss;
        oss << "[HaplotypeGrouping_IDRange] Chromosome id " << chromosome_id << " out of range.";
        throw runtime_error(oss.str().c_str());
    }

    --it;
    
    unsigned int id_start = *it;
    unsigned int group_index = it - id_starts_.begin();
    unsigned int id_count = id_counts_[group_index];

    if (chromosome_id >= id_start + id_count)
    {
        ostringstream oss;
        oss << "[HaplotypeGrouping_IDRange] Chromosome id " << chromosome_id << " out of range.";
        throw runtime_error(oss.str().c_str());
    }

    return group_index;
}


Parameters HaplotypeGrouping_IDRange::parameters() const
{
    Parameters parameters;

    for (vector<unsigned int>::const_iterator id_start=id_starts_.begin(), id_count=id_counts_.begin();
         id_start!=id_starts_.end(); ++id_start, ++id_count)
    {
        ostringstream value;
        value << *id_start << " " << *id_count;
        parameters.insert_name_value("start:count", value.str());
    }

    return parameters;
}


void HaplotypeGrouping_IDRange::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("start:count");
    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        istringstream iss(*it);
        unsigned int id_start = 0, id_count = 0;
        iss >> id_start >> id_count;
        id_starts_.push_back(id_start);
        id_counts_.push_back(id_count);
    }
}


//
// HaplotypeGrouping_Uniform
//


HaplotypeGrouping_Uniform::HaplotypeGrouping_Uniform(const string& id)
:   HaplotypeGrouping(id), ids_per_group_(0), id_count_(0)
{}


size_t HaplotypeGrouping_Uniform::group_count() const
{
    if (id_count_ == 0 || ids_per_group_ == 0)
        throw runtime_error("[HaplotypeGrouping_Uniform] Not initialized properly.");

    return id_count_ / ids_per_group_;
}


size_t HaplotypeGrouping_Uniform::group(unsigned int chromosome_id) const
{
    if (id_count_ == 0 || ids_per_group_ == 0)
        throw runtime_error("[HaplotypeGrouping_Uniform] Not initialized.");

    return (chromosome_id - id_offset_) / ids_per_group_;
}


Parameters HaplotypeGrouping_Uniform::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("ids_per_group", ids_per_group_);
    return parameters;
}


void HaplotypeGrouping_Uniform::configure(const Parameters& parameters, const Registry& registry)
{
    ids_per_group_ = parameters.value<size_t>("ids_per_group");
}


void HaplotypeGrouping_Uniform::initialize(const SimulatorConfig& config)
{
    if (!config.population_config_generator.get())
        throw runtime_error("[HaplotypeGrouping_Uniform] Null PopulationConfigGenerator.");

    Population::Configs popconfigs = 
        config.population_config_generator->population_configs(0, PopulationDataPtrs());

    if (popconfigs.size() != 1)
        throw runtime_error("[HaplotypeGrouping_Uniform] Invalid multiple populations.");

    id_offset_ = popconfigs[0].id_offset;
    id_count_ = 2 * popconfigs[0].population_size;
}


//
// Reporter_HaplotypeFrequencies
//


Reporter_HaplotypeFrequencies::Reporter_HaplotypeFrequencies(
    const string& id, 
    HaplotypeGroupingPtr haplotype_grouping,
    size_t chromosome_step,
    size_t update_step)
:   Configurable(id), 
    haplotype_grouping_(haplotype_grouping),
    chromosome_step_(chromosome_step),
    update_step_(update_step)
{}


void Reporter_HaplotypeFrequencies::update(size_t generation_index,
                                           const PopulationPtrs& populations,
                                           const PopulationDataPtrs& population_datas,
                                           bool is_final_generation)
{
    if (output_directory_.string().empty())
        throw runtime_error("[Reporter_HaplotypeFrequencies] No output directory specified.");

    if (!is_final_generation && (update_step_ == 0 || generation_index%update_step_ != 0)) return;

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_HaplotypeFrequencies] Population data size mismatch.");

    if (populations.empty()) return;

    for (size_t population_index=0; population_index<populations.size(); ++population_index)
    {
        const Population& population = *populations[population_index];

        for (size_t chromosome_pair_index=0; chromosome_pair_index!=chromosome_lengths_.size(); ++chromosome_pair_index)
        {
            ostringstream filename;
            filename << "haplotype_frequencies"
                    << "_chr" << chromosome_pair_index + 1;
            if (is_final_generation) 
                filename << "_final";
            else
                filename << "_gen" << generation_index;
            filename << "_pop" << population_index + 1
                     << ".txt";

            bfs::ofstream os(output_directory_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_HaplotypeFrequencies] Unable to open file " + filename.str()).c_str());

            write_haplotype_frequencies(os, population, chromosome_pair_index,
                                        chromosome_lengths_[chromosome_pair_index]);
        }
    }
}


Parameters Reporter_HaplotypeFrequencies::parameters() const
{
    Parameters parameters;
    if (haplotype_grouping_.get())
        parameters.insert_name_value("haplotype_grouping", haplotype_grouping_->object_id());
    parameters.insert_name_value("chromosome_step", chromosome_step_);
    parameters.insert_name_value("update_step", update_step_);
    return parameters;
}


void Reporter_HaplotypeFrequencies::configure(const Parameters& parameters, const Registry& registry)
{
    haplotype_grouping_ = registry.get<HaplotypeGrouping>(parameters.value<string>("haplotype_grouping"));
    chromosome_step_ = parameters.value<size_t>("chromosome_step");
    update_step_ = parameters.value<size_t>("update_step", 0);
}


void Reporter_HaplotypeFrequencies::initialize(const SimulatorConfig& config)
{
    Reporter::initialize(config);

    if (!haplotype_grouping_.get())
        throw runtime_error("[Reporter_HaplotypeFrequencies] Null HaplotypeGrouping.");

    haplotype_grouping_->initialize(config);

    chromosome_lengths_ = config.population_config_generator->chromosome_lengths();
}


void Reporter_HaplotypeFrequencies::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    if (haplotype_grouping_.get())
        haplotype_grouping_->write_configuration(os, ids_written);
}


void Reporter_HaplotypeFrequencies::write_haplotype_frequencies(std::ostream& os, 
    const Population& population, size_t chromosome_pair_index, unsigned int chromosome_length) const
{
    // header

    os << "# position group1 [group2 ...]\n";

    // one line per position

    for (size_t position=0; position<chromosome_length; position+=chromosome_step_)
    {
        vector<double> counts(haplotype_grouping_->group_count());
        double count_total = 0;

        ChromosomePairRangeIterator it = population.begin();
        ChromosomePairRangeIterator end = population.end();
        for (; it!=end; ++it)
        {
            if (it->size() == 0)
                throw runtime_error("[Reporter_HaplotypeFrequencies] I am insane.");

            const ChromosomePair& p = *(it->begin() + chromosome_pair_index);

            unsigned int id1 = p.first.find_haplotype_chunk(position)->id;
            unsigned int id2 = p.second.find_haplotype_chunk(position)->id;

            ++counts[haplotype_grouping_->group(id1)];
            ++counts[haplotype_grouping_->group(id2)];
            count_total += 2;
        }

        // write the haplotype frequencies
        
        os << position << " ";

        transform(counts.begin(), counts.end(), 
                  ostream_iterator<double>(os, " "),
                  bind2nd(divides<double>(), count_total));

        os << endl;            
    }
}


//
// Reporter_Regions
//


string Reporter_Regions::Region::configuration() const
{
    ostringstream result;
    result << locus.object_id() << " " << length;
    return result.str();
}


Reporter_Regions::Reporter_Regions(const std::string& id, 
                                   const Regions& regions)
:   Configurable(id), 
    regions_(regions),
    ms_mapping_begin_(0), ms_mapping_end_(0)
{}


void Reporter_Regions::update(size_t generation_index,
                              const PopulationPtrs& populations,
                              const PopulationDataPtrs& population_datas,
                              bool is_final_generation)
{
    if (!is_final_generation) return;

    if (output_directory_.string().empty())
        throw runtime_error("[Reporter_Regions] No output directory specified.");

    if (populations.size() != population_datas.size())
        throw runtime_error("[Reporter_Regions] Population data size mismatch.");

    if (populations.empty()) return;

    PopulationDataPtrs::const_iterator population_data = population_datas.begin();

    Loci loci_all;
    
    for (GenotypeMap::const_iterator it=(*population_data)->genotypes->begin(); 
         it!=(*population_data)->genotypes->end(); ++it)
        loci_all.insert(it->first);

    Loci loci_regions;

    for (Regions::const_iterator it=regions_.begin(); it!=regions_.end(); ++it)
    {
        Loci::const_iterator begin = loci_all.lower_bound(it->locus);
        Loci::const_iterator end = loci_all.lower_bound(Locus("", 
                                                              it->locus.chromosome_pair_index, 
                                                              it->locus.position + it->length));
        loci_regions.insert(begin, end);
    }

    size_t population_size_total = 0;
    map<Locus,size_t> variant_count_totals;

    for (size_t population_index=0; population_index<populations.size(); ++population_index, ++population_data)
    {
        ostringstream filename;
        filename << "regions" // TODO: change to object id for filename base
                << "_pop" << population_index + 1
                << ".txt";

        bfs::ofstream os(output_directory_ / filename.str());
        if (!os)
            throw runtime_error(("[Reporter_Regions] Unable to open file " + filename.str()).c_str());

        os << "# loci: ";
        copy(loci_regions.begin(), loci_regions.end(), ostream_iterator<Locus>(os, " "));
        os << endl;

        size_t population_size = populations[population_index]->population_size();
        population_size_total += population_size;
   
        // one line per individual

        for (size_t i=0; i<population_size; ++i)
        {
            os << i << ":";
            for (Loci::const_iterator locus=loci_regions.begin(); locus!=loci_regions.end(); ++locus)
                os << int(genotype_sum((*(*population_data)->genotypes->get(*locus))[i]));
            os << endl;
        }

        os.close();

        // sfs

        vector<size_t> sfs(2*population_size + 1);

        for (Loci::const_iterator locus=loci_regions.begin(); locus!=loci_regions.end(); ++locus)
        {
            const GenotypeData& genotypes = *(*population_data)->genotypes->get(*locus);
            size_t variant_count = 0;
            for (GenotypeData::const_iterator g=genotypes.begin(); g!=genotypes.end(); ++g)
                variant_count += size_t(genotype_sum(*g));
            if (variant_count > sfs.size())
                throw runtime_error("[Reporter_Regions] Only implemented for 0/1 variants");

            ++sfs[variant_count];
            variant_count_totals[*locus] += variant_count;
        }

        ostringstream filename_sfs;
        filename_sfs << "regions"
                     << "_pop" << population_index + 1
                     << ".sfs.txt";

        bfs::ofstream os_sfs(output_directory_ / filename_sfs.str());
        if (!os_sfs)
            throw runtime_error(("[Reporter_Regions] Unable to open file " + filename_sfs.str()).c_str());

        os_sfs << "# N = " << population_size << endl;
        copy(sfs.begin(), sfs.end(), ostream_iterator<size_t>(os_sfs, " "));
        os_sfs << endl;
        os_sfs.close();
    }

    // ms format output

    Loci loci_polymorphic;
    for (Loci::const_iterator locus=loci_regions.begin(); locus!=loci_regions.end(); ++locus)
    {
        unsigned int variant_count_total = variant_count_totals[*locus];
        if (variant_count_total > 0 && variant_count_total < population_size_total)
            loci_polymorphic.insert(*locus);
    }

    population_data = population_datas.begin();
    for (size_t population_index=0; population_index<populations.size(); ++population_index, ++population_data)
    {
        size_t population_size = populations[population_index]->population_size();

        ostringstream filename;
        filename << "regions" // TODO: change to object id for filename base
                << "_pop" << population_index + 1
                << ".ms";

        bfs::ofstream os(output_directory_ / filename.str());
        if (!os)
            throw runtime_error(("[Reporter_Regions] Unable to open file " + filename.str()).c_str());

        os << "ms " << population_size*2 << " 1\n";
        os << "0 0 0  # ms output generated by forqs [Reporter_Regions]\n";
        os << endl;
        
        os << "//\n";
        os << "segsites: " << loci_polymorphic.size() << endl;

        os << "positions: ";
        const double ms_mapping_length = ms_mapping_end_ - ms_mapping_begin_;
        for (Loci::const_iterator locus=loci_polymorphic.begin(); locus!=loci_polymorphic.end(); ++locus)
            os << (locus->position - ms_mapping_begin_)/ms_mapping_length << " ";
        os << endl;
   
        // two lines per individual

        for (size_t i=0; i<population_size; ++i)
        {
            // note: inefficient iteration

            for (Loci::const_iterator locus=loci_polymorphic.begin(); locus!=loci_polymorphic.end(); ++locus)
                os << int(genotype_first((*(*population_data)->genotypes->get(*locus))[i]));
            os << endl;

            for (Loci::const_iterator locus=loci_polymorphic.begin(); locus!=loci_polymorphic.end(); ++locus)
                os << int(genotype_second((*(*population_data)->genotypes->get(*locus))[i]));
            os << endl;
        }

        os.close();
    }
}
    

Parameters Reporter_Regions::parameters() const
{
    Parameters parameters;
    for (Regions::const_iterator it=regions_.begin(); it!=regions_.end(); ++it)
        parameters.insert_name_value("locus:length", it->configuration());
    parameters.insert_name_value("ms_mapping_begin", ms_mapping_begin_);
    parameters.insert_name_value("ms_mapping_end", ms_mapping_end_);
    return parameters;
}


void Reporter_Regions::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("locus:length");

    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        istringstream iss(*it);
        string id_locus;
        size_t length = 0;
        iss >> id_locus >> length;

        if (id_locus.empty() || !iss)
            throw runtime_error(("[Reporter_Regions] Error parsing locus:length: " + *it).c_str());

        regions_.push_back(Region(*registry.get<Locus>(id_locus), length));
    }

    // regions must be on same chromosome for ms-format mapping to make sense

    if (regions_.empty())
        throw runtime_error("[Reporter_Regions] No regions specified.");

    size_t chromosome_pair_index = regions_.front().locus.chromosome_pair_index;
    for (Regions::const_iterator it=regions_.begin()+1; it!=regions_.end(); ++it)
        if (it->locus.chromosome_pair_index != chromosome_pair_index)
            throw runtime_error("[Reporter_Regions] Regions must be on same chromosome.");

    // if mapping parameters not specified, map smallest enclosing interval [begin,end] -> [0,1]

    ms_mapping_begin_ = parameters.value<unsigned int>("ms_mapping_begin", numeric_limits<unsigned int>::max());
    ms_mapping_end_ = parameters.value<unsigned int>("ms_mapping_end", 0);

    if (ms_mapping_end_ == 0)
    {
        for (Regions::const_iterator it=regions_.begin(); it!=regions_.end(); ++it)
        {
            if (it->locus.position < ms_mapping_begin_)
                ms_mapping_begin_ = it->locus.position;

            if (it->locus.position + it->length > ms_mapping_end_)
                ms_mapping_end_ = it->locus.position + it->length;
        }
    }
}


void Reporter_Regions::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    for (Regions::const_iterator it=regions_.begin(); it!=regions_.end(); ++it)
        it->locus.write_configuration(os, ids_written);
}


//
// Reporter_DeterministicTrajectories
//


void Reporter_DeterministicTrajectories::update(size_t generation_number,
                    const PopulationPtrs& populations,
                    const PopulationDataPtrs& population_datas,
                    bool is_final_generation)
{
    if (!is_final_generation) return;

    bfs::ofstream os_p(output_directory_ / "deterministic_trajectories_allele_frequencies.txt");
    bfs::ofstream os_w(output_directory_ / "deterministic_trajectories_mean_fitnesses.txt");

    double p = initial_allele_frequency_;
    
    for (size_t generation=0; generation<=generation_count_; ++generation)
    {
        double q = 1 - p;
        double w0_marginal = q*w_[0] + p*w_[1];
        double w1_marginal = q*w_[1] + p*w_[2];
        double w_mean = q*w0_marginal + p*w1_marginal;
        double dp = p*q*(w1_marginal-w0_marginal)/w_mean;
        //double dw = dp*(w1_marginal-w0_marginal) + dp*dp*(w2-2*w1+w0);

        os_p << p << endl;
        os_w << w_mean << endl;

        p = p + dp;
    }
}


Parameters Reporter_DeterministicTrajectories::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("initial_allele_frequency", initial_allele_frequency_);
    parameters.insert_name_value("w0", w_[0]);
    parameters.insert_name_value("w1", w_[1]);
    parameters.insert_name_value("w2", w_[2]);
    return parameters;
}


void Reporter_DeterministicTrajectories::configure(const Parameters& parameters, const Registry& registry)
{
    initial_allele_frequency_ = parameters.value<double>("initial_allele_frequency");

    if (parameters.count("additive_selection_coefficient"))
    {
        double s = parameters.value<double>("additive_selection_coefficient");
        w_[0] = 1;
        w_[1] = 1 + s;
        w_[2] = 1 + 2*s;
    }
    else
    {
        w_[0] = parameters.value<double>("w0");
        w_[1] = parameters.value<double>("w1");
        w_[2] = parameters.value<double>("w2");
    }
}


void Reporter_DeterministicTrajectories::initialize(const SimulatorConfig& config)
{
    generation_count_ = config.population_config_generator->generation_count();
    Reporter::initialize(config);
}


//
// Reporter_Variants
//


Reporter_Variants::Reporter_Variants(const string& id)
:   Configurable(id), update_step_(0)
{}


void Reporter_Variants::update(size_t generation_index,
                            const PopulationPtrs& populations,
                            const PopulationDataPtrs& population_datas,
                            bool is_final_generation)
{
    if (update_step_ && (generation_index % update_step_ == 0) || 
        generation_index == 0 ||
        is_final_generation)
    {
        if (populations.size() != population_datas.size())
            throw runtime_error("[Reporter_Variants] Population data size mismatch.");

        const size_t population_count = populations.size();
        const char* filestem = "variants";

        for (size_t population_index=0; population_index<population_count; ++population_index)
        {
            ostringstream filename;

            if (is_final_generation)
                filename << filestem << "_final_pop" << population_index + 1 << ".txt"; 
            else
                filename << filestem << "_gen" << generation_index << "_pop" << population_index + 1 << ".txt"; 

            bfs::ofstream os(output_directory_ / filename.str());
            if (!os)
                throw runtime_error(("[Reporter_Variants] Unable to open " + filename.str()).c_str());

            const PopulationData& data = *population_datas[population_index];
            const size_t population_size = data.population_size;

            // iterate through genotype data in parallel (output columns)

            vector<GenotypeData::const_iterator> its;

            for (vector<Locus>::const_iterator locus=loci_.begin(); locus!=loci_.end(); ++locus)
            {
                const GenotypeData& genotypes = *data.genotypes->get(*locus);
                if (genotypes.size() != population_size)
                    throw runtime_error("[Reporter_Variants] Data size mismatch.");
                its.push_back(genotypes.begin());
            }

            os << "#positions: ";
            copy(loci_.begin(), loci_.end(), ostream_iterator<Locus>(os, " "));
            os << endl;
            os << "segsites: " << loci_.size() << endl;
            
            for (size_t i=0; i<population_size; ++i)
            {
                // write first haplotype

                for (vector<GenotypeData::const_iterator>::iterator it=its.begin(); 
                     it!=its.end(); ++it)
                {
                    os << int(genotype_first(**it));
                }
                os << endl;

                // write second haplotype

                for (vector<GenotypeData::const_iterator>::iterator it=its.begin(); 
                     it!=its.end(); ++it)
                {
                    os << int(genotype_second(**it));
                    ++(*it);
                }
                os << endl;
            }

            os.close();
        }
    }
}


Loci Reporter_Variants::loci(size_t generation_index, bool is_final_generation) const
{
    Loci loci;
    copy(loci_.begin(), loci_.end(), inserter(loci, loci.begin()));
    return loci;
}


Parameters Reporter_Variants::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value_vector("loci", locus_ids_);
    parameters.insert_name_value("update_step", update_step_);
    return parameters;
}


void Reporter_Variants::configure(const Parameters& parameters, const Registry& registry)
{
    locus_ids_ = parameters.value_vector<string>("loci");
    update_step_ = parameters.value<size_t>("update_step", 0);

    for (vector<string>::const_iterator id=locus_ids_.begin(); id!=locus_ids_.end(); ++id)
    {
        LocusPtr locus = registry.get<Locus>(*id, std::nothrow);
        LocusListPtr locus_list = registry.get<LocusList>(*id, std::nothrow);
        QuantitativeTraitPtr qt = registry.get<QuantitativeTrait>(*id, std::nothrow);

        if (locus.get()) 
        {
            loci_.push_back(*locus);
            children_.push_back(dynamic_pointer_cast<Configurable>(locus));
        }
        else if (locus_list.get())
        {
            copy(locus_list->begin(), locus_list->end(), back_inserter(loci_));
            children_.push_back(dynamic_pointer_cast<Configurable>(locus_list));
        }
        else if (qt.get())
        {
            qts_specified_.push_back(qt); // wait until initialize() to add to loci_
            children_.push_back(dynamic_pointer_cast<Configurable>(qt));
        }
        else
            throw runtime_error("[VariantIndicator_File] id must be Locus or LocusList: " + *id);
    }
}


void Reporter_Variants::initialize(const SimulatorConfig& config)
{
    Reporter::initialize(config);
    for (QuantitativeTraitPtrs::const_iterator it=qts_specified_.begin(); it!=qts_specified_.end(); ++it)
        copy((*it)->loci().begin(), (*it)->loci().end(), inserter(loci_, loci_.end()));
}


void Reporter_Variants::write_child_configurations(ostream& os, set<string>& written) const
{
    for (ConfigurablePtrs::const_iterator it=children_.begin(); it!=children_.end(); ++it)
        (*it)->write_configuration(os, written);
}


