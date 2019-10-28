//
// VariantIndicatorImplementation.cpp
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


#include "VariantIndicatorImplementation.hpp"


using namespace std;


//
// VariantIndicator_Composite
//


unsigned int VariantIndicator_Composite::operator()(unsigned int chunk_id, const Locus& locus) const
{
    for (VariantIndicatorPtrs::const_iterator it=variant_indicators_.begin(); it!=variant_indicators_.end(); ++it)
    {
        unsigned int value = (**it)(chunk_id, locus);
        if (value != 0) return value;
    }

    return 0;
}


Parameters VariantIndicator_Composite::parameters() const 
{
    ostringstream ids;
    for (VariantIndicatorPtrs::const_iterator it=variant_indicators_.begin(); it!=variant_indicators_.end(); ++it)
        ids << (*it)->object_id() << " ";

    Parameters parameters;
    parameters.insert_name_value("variant_indicators", ids.str());
    return parameters;
}


void VariantIndicator_Composite::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> vi_ids = parameters.value_vector<string>("variant_indicators");
    for (vector<string>::const_iterator id=vi_ids.begin(); id!=vi_ids.end(); ++id)
        variant_indicators_.push_back(registry.get<VariantIndicator>(*id));
}


void VariantIndicator_Composite::write_child_configurations(ostream& os, set<string>& written) const
{
    for (VariantIndicatorPtrs::const_iterator it=variant_indicators_.begin(); it!=variant_indicators_.end(); ++it)
        (*it)->write_configuration(os, written);
}


//
// VariantIndicator_IDRange
//


unsigned int VariantIndicator_IDRange::operator()(unsigned int chunk_id, const Locus& locus) const
{
    pair<EntryMap::const_iterator,EntryMap::const_iterator> range = entries_.equal_range(locus);
    for (EntryMap::const_iterator it=range.first; it!=range.second; ++it)
    {
        const Entry& entry = it->second;
        if (chunk_id >= entry.id_start && 
            chunk_id < entry.id_start + entry.id_count &&
            ((chunk_id - entry.id_start) % entry.id_step) == 0)
            return entry.value;
    }
    return 0;
}


Parameters VariantIndicator_IDRange::parameters() const 
{
    Parameters parameters;

    for (EntryMap::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
    {
        ostringstream oss;
        oss << it->first.object_id() << " "
            << it->second.id_start << " "
            << it->second.id_count << " "
            << it->second.id_step << " "
            << it->second.value;
        parameters.insert_name_value("locus:start:count:step:value", oss.str());
    }

    return parameters;
}


void VariantIndicator_IDRange::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("locus:start:count:step:value");
    
    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        string id_locus;
        unsigned int id_start = 0, id_count = 0, id_step = 0, value = 0;

        istringstream iss(*it);
        iss >> id_locus >> id_start >> id_count >> id_step >> value;

        if (id_step == 0)
            throw runtime_error("[VariantIndicator_IDRange] Error: id_step == 0"); 
    
        const Locus& locus = *registry.get<Locus>(id_locus);

        entries_.insert(make_pair(locus, Entry(id_start, id_count, id_step, value)));
    }
}


void VariantIndicator_IDRange::write_child_configurations(ostream& os, set<string>& written) const
{
    for (EntryMap::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
            it->first.write_configuration(os, written);
}


//
// VariantIndicator_IDSet
//


unsigned int VariantIndicator_IDSet::operator()(unsigned int chunk_id, const Locus& locus) const
{
    if (!entries_.count(locus)) return 0;
    const Entry& entry = entries_.at(locus);
    return entry.ids.count(chunk_id) ? entry.value: 0;
}


void VariantIndicator_IDSet::write_file(const std::string& filename) const
{
    ofstream os(filename.c_str());

    for (EntryMap::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
    {
        os << "locus:value:ids = " 
           << it->first.object_id() << " "
           << it->second.value << " ";
        copy(it->second.ids.begin(), it->second.ids.end(), ostream_iterator<unsigned int>(os, " "));
        os << endl;
    }
}


Parameters VariantIndicator_IDSet::parameters() const 
{
    Parameters parameters;

    for (EntryMap::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
    {
        ostringstream oss;
        oss << it->first.object_id() << " "
            << it->second.value << " ";
        copy(it->second.ids.begin(), it->second.ids.end(), ostream_iterator<double>(oss, " "));
        parameters.insert_name_value("locus:value:ids", oss.str());
    }

    return parameters;
}


void VariantIndicator_IDSet::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("locus:value:ids");
    
    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        string id_locus;
        unsigned int value = 0;
        set<unsigned int> ids;

        istringstream iss(*it);
        iss >> id_locus >> value;
        copy(istream_iterator<unsigned int>(iss), istream_iterator<unsigned int>(), inserter(ids, ids.end()));

        const Locus& locus = *registry.get<Locus>(id_locus);

        entries_.insert(make_pair(locus, Entry(value, ids)));
    }
}


void VariantIndicator_IDSet::write_child_configurations(ostream& os, set<string>& written) const
{
    for (EntryMap::const_iterator it=entries_.begin(); it!=entries_.end(); ++it)
        it->first.write_configuration(os, written);
}


//
// VariantIndicator_Random
//


//
// implementation note:
//      (Info::population_index == -1ul) == "all populations"
//


Parameters VariantIndicator_Random::parameters() const 
{
    Parameters parameters;

    for (LocusListInfos::const_iterator it=locus_list_infos_.begin(); it!=locus_list_infos_.end(); ++it)
    {
        ostringstream oss;
        oss << it->locus_list->object_id() << " ";

        if (it->population_index == -1ul)
            oss << "* ";
        else
            oss << it->population_index + 1 << " "; // 1-based

        if (it->distribution.get())
        {
            oss << it->distribution->object_id();
            parameters.insert_name_value("locus_list:population:frequency_distribution", oss.str());
        }
        else
        {
            copy(it->frequencies.begin(), it->frequencies.end(), ostream_iterator<double>(oss, " "));
            parameters.insert_name_value("locus_list:population:frequencies", oss.str());
        }
    }

    return parameters;
}


void VariantIndicator_Random::configure(const Parameters& parameters, const Registry& registry)
{
    vector<string> configurations = parameters.values<string>("locus_list:population:frequencies");
    
    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        // locus_list:population:frequencies = id population freq1 freq2 ... freqN

        string id_locus_list, population_string;
        vector<double> frequencies;

        istringstream iss(*it);
        iss >> id_locus_list >> population_string;
        copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(frequencies));

        if (frequencies.empty())
            throw runtime_error("[VariantIndicator_Random] Empty frequencies.");

        LocusListPtr locus_list = registry.get<LocusList>(id_locus_list);

        if (population_string == "*") // all populations
        {
            locus_list_infos_.push_back(LocusListInfo(-1ul, locus_list, frequencies));
        }
        else
        {
            unsigned int population = atoi(population_string.c_str()); // 1-based
            if (population == 0) throw runtime_error("[VariantIndicator_Random] Invalid population 0.");
            locus_list_infos_.push_back(LocusListInfo(population-1, locus_list, frequencies));
        }
    }

    configurations = parameters.values<string>("locus_list:population:frequency_distribution");

    for (vector<string>::const_iterator it=configurations.begin(); it!=configurations.end(); ++it)
    {
        // locus_list:population:frequency_distribution = id population id_dist

        string id_locus_list, population_string, id_distribution;

        istringstream iss(*it);
        iss >> id_locus_list >> population_string >> id_distribution;

        LocusListPtr locus_list = registry.get<LocusList>(id_locus_list);
        Random::DistributionPtr distribution = registry.get<Random::Distribution>(id_distribution);

        if (population_string == "*") // all populations
        {
            locus_list_infos_.push_back(LocusListInfo(-1ul, locus_list, distribution));
        }
        else
        {
            unsigned int population = atoi(population_string.c_str()); // 1-based
            if (population == 0) throw runtime_error("[VariantIndicator_Random] Invalid population 0.");
            locus_list_infos_.push_back(LocusListInfo(population-1, locus_list, distribution));
        }
    }
}


void VariantIndicator_Random::initialize(const SimulatorConfig& simconfig)
{
    if (!simconfig.population_config_generator.get())
        throw runtime_error("[VariantIndicator_Random] Null population_config_generator.");

    Population::Configs popconfigs = 
        simconfig.population_config_generator->population_configs(0, PopulationDataPtrs());

    for (Population::Configs::const_iterator popconfig=popconfigs.begin(); popconfig!=popconfigs.end(); ++popconfig)
    {
        const size_t population_index = popconfig - popconfigs.begin();
        const unsigned int id_start = popconfig->id_offset;
        const unsigned int id_count = 2 * popconfig->population_size;
        const unsigned int value = 1;

        for (LocusListInfos::iterator info=locus_list_infos_.begin(); info!=locus_list_infos_.end(); ++info)
        {
            if (info->population_index != -1ul && info->population_index >= popconfigs.size())
                throw runtime_error("[VariantIndicator_Random] Bad population index.");

            if (info->population_index != -1ul && info->population_index != population_index)
                continue;

            if (info->frequencies.empty() && info->distribution.get())
            {
                for (size_t i=0; i<info->locus_list->size(); ++i)
                    info->frequencies.push_back(info->distribution->random_value());
            }

            if (info->frequencies.size() != info->locus_list->size())
            {
                cout << "frequencies.size(): " << info->frequencies.size() << endl;
                cout << "locus_list->size(): " << info->locus_list->size() << endl;
                throw runtime_error("[VariantIndicator_Random] Frequency count does not match LocusList size.");
            }

            vector<double>::const_iterator f = info->frequencies.begin();
            for (LocusList::const_iterator locus=info->locus_list->begin(); 
                 locus!=info->locus_list->end(); ++locus, ++f)
            {
                vector<size_t> indices = 
                    Random::random_indices_without_replacement(id_count, size_t(*f * id_count));

                set<unsigned int> ids;
                for (vector<size_t>::const_iterator index=indices.begin(); index!=indices.end(); ++index)
                    ids.insert(id_start + *index);
                
                if (entries_.count(*locus) == 0)
                    entries_.insert(make_pair(*locus, Entry(value)));
                copy(ids.begin(), ids.end(), 
                     inserter(entries_.at(*locus).ids, entries_.at(*locus).ids.begin()));
            }
        }
    }
}


void VariantIndicator_Random::write_child_configurations(ostream& os, set<string>& written) const
{
    for (LocusListInfos::const_iterator it=locus_list_infos_.begin(); it!=locus_list_infos_.end(); ++it)
    {
        it->locus_list->write_configuration(os, written);
        if (it->distribution.get())
            it->distribution->write_configuration(os, written);
    }
}


//
// VariantIndicator_File
//


unsigned int VariantIndicator_File::operator()(unsigned int chunk_id, const Locus& locus) const
{
    if (!ms_.get())
        throw runtime_error("[VariantIndicator_File] Null ms file pointer.");
    
    if (chunk_id >= ms_->sequences.size())
    {
        ostringstream oss;
        oss << "[VariantIndicator_File] Haplotype id " << chunk_id << " out of range.";
        throw runtime_error(oss.str().c_str());
    }

    if (locus_index_map_.count(locus))
    {
        size_t locus_index = locus_index_map_.at(locus);
        return ms_->sequences[chunk_id][locus_index];
    }

    return 0;
}


Parameters VariantIndicator_File::parameters() const 
{
    Parameters parameters;
    parameters.insert_name_value("msfile", msfile_);
    parameters.insert_name_value_vector("loci", locus_ids_);
    return parameters;
}


void VariantIndicator_File::configure(const Parameters& parameters, const Registry& registry)
{
    msfile_ = parameters.value<string>("msfile");
    ms_ = MSFormatPtr(new MSFormat(msfile_));
    locus_ids_ = parameters.value_vector<string>("loci");

    if (locus_ids_.size() > ms_->segsites())
        throw runtime_error("[VariantIndicator_File] More loci specified than segsites in ms file.");

    for (vector<string>::const_iterator id=locus_ids_.begin(); id!=locus_ids_.end(); ++id)
    {
        LocusPtr locus = registry.get<Locus>(*id, std::nothrow);
        LocusListPtr locus_list = registry.get<LocusList>(*id, std::nothrow);

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
        else
            throw runtime_error("[VariantIndicator_File] id must be Locus or LocusList:" + *id);
    }

    for (vector<Locus>::const_iterator it=loci_.begin(); it!=loci_.end(); ++it)
        locus_index_map_[*it] = it-loci_.begin();

    // fix up ms sequences (char -> int)
    for (vector<string>::iterator it=ms_->sequences.begin(); it!=ms_->sequences.end(); ++it)
        for (string::iterator jt=it->begin(); jt!=it->end(); ++jt)
            *jt = char(boost::lexical_cast<int>(*jt));
}


void VariantIndicator_File::write_child_configurations(ostream& os, set<string>& written) const
{
    for (ConfigurablePtrs::const_iterator it=children_.begin(); it!=children_.end(); ++it)
        (*it)->write_configuration(os, written);
}


//
// VariantIndicator_SingleLocusHardyWeinberg
//


VariantIndicator_SingleLocusHardyWeinberg::VariantIndicator_SingleLocusHardyWeinberg(const string& id, 
                                                                                     Locus locus,
                                                                                     double allele_frequency)
:   Configurable(id), // ?
    VariantIndicator_IDRange(id),
    locus_(locus),
    allele_frequency_(allele_frequency)
{}


Parameters VariantIndicator_SingleLocusHardyWeinberg::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("locus", locus_.object_id());
    parameters.insert_name_value("allele_frequency", allele_frequency_);
    return parameters;
}


void VariantIndicator_SingleLocusHardyWeinberg::configure(const Parameters& parameters, const Registry& registry)
{
    locus_ = *registry.get<Locus>(parameters.value<string>("locus"));
    allele_frequency_ = parameters.value<double>("allele_frequency");
}


void VariantIndicator_SingleLocusHardyWeinberg::initialize(const SimulatorConfig& simconfig)
{
    if (!simconfig.population_config_generator.get())
        throw runtime_error("[VariantIndicator_SingleLocusHardyWeinberg] Null population_config_generator.");

    Population::Configs popconfigs = 
        simconfig.population_config_generator->population_configs(0, PopulationDataPtrs());

    const double p = allele_frequency_;
    const double q = 1-p;
    const unsigned int value = 1;

    for (Population::Configs::const_iterator popconfig=popconfigs.begin(); popconfig!=popconfigs.end(); ++popconfig)
    {
        // individuals: (0, 1, 2,                 ...                , 2N-1)
        // genotypes:   (2, 2, 2, ... , 2,  1, 1,  ... , 1,  0, 0,  ... , 0)
        //                       N*p^2             N(2pq)           N*q^2
        // note: 2 consecutive ids assigned per individual

        const unsigned int id_start = popconfig->id_offset;
        const unsigned int id_count = 2 * popconfig->population_size;

        const unsigned int count2 = (unsigned int)round(id_count*p*p);
        const unsigned int count1 = (unsigned int)round(id_count*(2*p*q));

        entries_.insert(make_pair(locus_, Entry(id_start, count2, 1, value))); // homozygote derived
        entries_.insert(make_pair(locus_, Entry(id_start+count2, count1, 2, value))); // heterozygote
    }
}


void VariantIndicator_SingleLocusHardyWeinberg::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    locus_.write_configuration(os, ids_written);
}


// 
// VariantIndicator_TwoLocusLD 
//


VariantIndicator_TwoLocusLD::VariantIndicator_TwoLocusLD(const string& id, 
                            size_t population_size,
                            Locus locus_1,
                            double allele_frequency_1,
                            Locus locus_2,
                            double allele_frequency_2,
                            double D,
                            unsigned int id_offset_step)
:   Configurable(id), 
    population_size_(population_size),
    locus_1_(locus_1),
    allele_frequency_1_(allele_frequency_1),
    locus_2_(locus_2),
    allele_frequency_2_(allele_frequency_2),
    D_(D),
    id_offset_step_(id_offset_step)
{ 
    initialize_internal();
}


unsigned int VariantIndicator_TwoLocusLD::operator()(unsigned int chunk_id, const Locus& locus) const
{
    if (locus != locus_1_ && locus != locus_2_) return 0;

    unsigned int id = id_offset_step_ ? chunk_id%id_offset_step_ : chunk_id;

    if (id < max_00_) return 0;
    else if (id < max_01_) return locus == locus_2_;
    else if (id < max_10_) return locus == locus_1_;
    else if (id < max_11_) return 1;

    return 0; 
}


Parameters VariantIndicator_TwoLocusLD::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("population_size", population_size_);
    parameters.insert_name_value("locus_1", locus_1_.object_id());
    parameters.insert_name_value("allele_frequency_1", allele_frequency_1_);
    parameters.insert_name_value("locus_2", locus_2_.object_id());
    parameters.insert_name_value("allele_frequency_2", allele_frequency_2_);
    parameters.insert_name_value("D", D_);
    parameters.insert_name_value("id_offset_step", id_offset_step_);
    return parameters;
}


void VariantIndicator_TwoLocusLD::configure(const Parameters& parameters, const Registry& registry)
{
    population_size_ = parameters.value<size_t>("population_size");
    locus_1_ = *registry.get<Locus>(parameters.value<string>("locus_1"));
    allele_frequency_1_ = parameters.value<double>("allele_frequency_1");
    locus_2_ = *registry.get<Locus>(parameters.value<string>("locus_2"));
    allele_frequency_2_ = parameters.value<double>("allele_frequency_2");
    D_ = parameters.value<double>("D");
    id_offset_step_ = parameters.value<size_t>("id_offset_step", 0);

    initialize_internal();
}


void VariantIndicator_TwoLocusLD::write_child_configurations(ostream& os, set<string>& ids_written) const
{
    locus_1_.write_configuration(os, ids_written);
    locus_2_.write_configuration(os, ids_written);
}


void VariantIndicator_TwoLocusLD::initialize_internal()
{
    const size_t N = population_size_;
    const double f1 = allele_frequency_1_;
    const double f2 = allele_frequency_2_;

    const double p00 = f1*f2 + D_;
    const double p01 = (1-f1)*f2 - D_;
    const double p10 = f1*(1-f2) - D_;
    const double p11 = (1-f1)*(1-f2) + D_;

    // note: 2 consecutive ids assigned per individual

    max_00_ = (size_t)round(2 * N * p00);
    max_01_ = max_00_ + (size_t)round(2 * N * p01);
    max_10_ = max_01_ + (size_t)round(2 * N * p10);
    max_11_ = max_10_ + (size_t)round(2 * N * p11);
}


//
// VariantIndicator_Mutable
//


VariantIndicator_Mutable::VariantIndicator_Mutable(const string& id, 
                                                   unsigned int unused_id_start, 
                                                   VariantIndicatorPtr vi,
                                                   const string& output_directory)
:   Configurable(id), unused_id_start_(unused_id_start), unused_id_current_(unused_id_start), 
    vi_(vi), outdir_(output_directory)
{
    const bool debug = false;
    if (debug && !output_directory.empty())
    {
        string filename = "debug_vi_mutable.txt";
        os_debug_.open(outdir_ / filename);
        if (!os_debug_)
            throw runtime_error(("[VariantIndicator_Mutable] Unable to open file " + filename).c_str());
    }
} 


unsigned int VariantIndicator_Mutable::operator()(unsigned int chunk_id, const Locus& locus) const
{   
    const IDValueMap empty;
    const IDValueMap& id_value_map = id_value_maps_.count(locus) ? id_value_maps_.at(locus) : empty;

    unsigned int parent_id = chunk_id;

    while (id_ancestry_.count(parent_id))
    {
        if (id_value_map.count(parent_id))
            return id_value_map.at(parent_id);  // return mutant value

        parent_id = id_ancestry_.at(parent_id); // walk back along chunk ancestry tree
    }

    return vi_.get() ? (*vi_)(parent_id, locus) : 0; // default to internal VariantIndicator
}


unsigned int VariantIndicator_Mutable::mutate(unsigned int old_chunk_id, const Locus& locus, unsigned int value)
{
    unsigned int new_chunk_id = unused_id_current_++;
    id_ancestry_[new_chunk_id] = old_chunk_id;
    id_value_maps_[locus][new_chunk_id] = value;
    return new_chunk_id;
}


Parameters VariantIndicator_Mutable::parameters() const
{
    Parameters parameters;
    parameters.insert_name_value("unused_id_start", unused_id_start_);
    if (vi_.get()) parameters.insert_name_value("variant_indicator", vi_->object_id());
    return parameters;
}


void VariantIndicator_Mutable::configure(const Parameters& parameters, const Registry& registry)
{
    unused_id_current_ = unused_id_start_ = parameters.value<unsigned int>("unused_id_start");

    if (parameters.count("variant_indicator"))
        vi_ = registry.get<VariantIndicator>(parameters.value<string>("variant_indicator"));
}


void VariantIndicator_Mutable::initialize(const SimulatorConfig& simconfig)
{
    vi_->initialize(simconfig);
}


void VariantIndicator_Mutable::write_child_configurations(ostream& os, set<string>& written) const
{
    if (vi_.get()) vi_->write_configuration(os, written);
}


void VariantIndicator_Mutable::report(std::ostream& os) const
{
    os << "[VariantIndicator_Mutable]\n";
    os << "unused_id_current_: " << unused_id_current_ << endl;

    os << "id_ancestry_:\n";
    for (IDAncestry::const_iterator it=id_ancestry_.begin(); it!=id_ancestry_.end(); ++it)
        os << "  " << it->first << " -> " << it->second << endl;

    os << "id_value_maps_\n";
    for (IDValueMaps::const_iterator it=id_value_maps_.begin(); it!=id_value_maps_.end(); ++it)
    {
        os << "  locus " << it->first << " values: ";
        for (IDValueMap::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
            os << jt->first << ":" << jt->second << " ";
        os << endl;
    }
}


void VariantIndicator_Mutable::update(size_t generation_index,
                                      const PopulationPtrs& populations,
                                      const PopulationDataPtrs& population_datas,
                                      bool is_final_generation)
{
    if (os_debug_.is_open())
    {
        os_debug_ << "generation " << generation_index << endl;
        report(os_debug_);
        os_debug_ << "\n\n";
    }

    if (is_final_generation)
    {
        bfs::ofstream os(outdir_ / "forqs.id_ancestry_map.txt");
        os << "# forqs id ancestry map [VariantIndicator_Mutable]\n";
        for (IDAncestry::const_iterator it=id_ancestry_.begin(); it!=id_ancestry_.end(); ++it)
            os << it->first << " " << it->second << endl;
        os.close();
    }
}


Loci VariantIndicator_Mutable::loci(size_t generation_index, 
                                    bool is_final_generation) const
{
    Loci loci;

    if (is_final_generation)
        for (IDValueMaps::const_iterator it=id_value_maps_.begin(); it!=id_value_maps_.end(); ++it)
            loci.insert(it->first);

    return loci;
}


