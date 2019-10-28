//
// forqs_focal_subset.cpp
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


#include "Population_Organisms.hpp"
#include "VariantIndicatorImplementation.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <numeric>


using namespace std;


struct Config
{
    string filename_full_population;
    string filename_full_ms;
    string filename_original_snps;
    string filename_subset_population;
    string filename_subset_ms;
    size_t population_size;
    LocusPtr focal_locus;
    double focal_allele_frequency;

    Config(const string& filename = "");

    void validate() const;
};


ostream& operator<<(ostream& os, const Config& config)
{
    os << "filename_full_population " << config.filename_full_population << endl
       << "filename_full_ms " << config.filename_full_ms << endl
       << "filename_original_snps " << config.filename_original_snps << endl
       << "filename_subset_population " << config.filename_subset_population << endl
       << "filename_subset_ms " << config.filename_subset_ms << endl
       << "population_size " << config.population_size << endl
       << "focal_locus " << config.focal_locus->chromosome_pair_index << " " << config.focal_locus->position << endl
       << "focal_allele_frequency " << config.focal_allele_frequency << endl;
    return os;
}


istream& operator>>(istream& is, Config& config)
{
    while(is)
    {
        string buffer;
        getline(is, buffer);
        if (!is) break;
        if (buffer.empty() || buffer[0] == '#') continue;
        istringstream iss(buffer);
        string parameter;
        iss >> parameter;

        if (parameter == "filename_full_population")
            iss >> config.filename_full_population;
        else if (parameter == "filename_full_ms")
            iss >> config.filename_full_ms;
        else if (parameter == "filename_original_snps")
            iss >> config.filename_original_snps;
        else if (parameter == "filename_subset_population")
            iss >> config.filename_subset_population;
        else if (parameter == "filename_subset_ms")
            iss >> config.filename_subset_ms;
        else if (parameter == "population_size")
            iss >> config.population_size;
        else if (parameter == "focal_locus")
            iss >> config.focal_locus->chromosome_pair_index >> config.focal_locus->position;
        else if (parameter == "focal_allele_frequency")
            iss >> config.focal_allele_frequency;
        else
            cerr << "Warning: unused parameter:\n" << buffer << endl;
    }
        
    return is;
}


Config::Config(const string& filename)
:   population_size(0), focal_locus(new Locus("focal_locus")), focal_allele_frequency(0)
{
    if (!filename.empty())
    {
        ifstream is(filename.c_str());
        is >> *this;
    }
}


void Config::validate() const
{
    if (filename_full_population.empty()) throw runtime_error("empty filename_full_population");
    if (filename_original_snps.empty()) throw runtime_error("empty filename_original_snps");
    if (filename_subset_population.empty()) throw runtime_error("empty filename_subset_population");
}


void create_subset_population(const Config& config)
{
    Population_Organisms population;
    ifstream is(config.filename_full_population.c_str());
    if (!is)
        throw runtime_error("Unable to open file " + config.filename_full_population);
    is >> population;
    is.close();

    cout << "full population size: " << population.population_size() << endl;

    VariantIndicatorPtr vi;

    if (config.focal_locus->chromosome_pair_index == 0 && config.focal_locus->position == 0)
    {
        vi = VariantIndicatorPtr(new VariantIndicator_Trivial("vi"));
    }
    else
    {
        vi = VariantIndicatorPtr(new VariantIndicator_File("vi_file"));

        Parameters parameters;
        parameters.insert_name_value("msfile", config.filename_original_snps);
        parameters.insert_name_value("loci", config.focal_locus->object_id());

        Configurable::Registry registry;
        registry[config.focal_locus->object_id()] = config.focal_locus;

        vi->configure(parameters, registry);
    }

    Genotyper genotyper;
    Loci loci;
    loci.insert(*config.focal_locus);
    GenotypeMap genotype_map;
    genotyper.genotype(loci, population, *vi, genotype_map);
    const GenotypeData& genotypes = *genotype_map.get(*config.focal_locus);

    vector< set<size_t> > genotype_bins(3);

    for (size_t i=0; i<genotypes.size(); ++i)
    {
        size_t genotype = size_t(genotype_sum(genotypes[i]));
        genotype_bins.at(genotype).insert(i);
    }

    vector<size_t> target_genotype_counts(3);
    const double p = config.focal_allele_frequency;
    const double q = 1-p;
    target_genotype_counts[0] = size_t(round(q*q*config.population_size));
    target_genotype_counts[1] = size_t(round(2*p*q*config.population_size));
    target_genotype_counts[2] = size_t(round(p*p*config.population_size));
    size_t count_diff = config.population_size - 
        accumulate(target_genotype_counts.begin(), target_genotype_counts.end(), 0);
    if (count_diff != 0) cerr << "Warning: Hardy-Weinberg not exact.";
    target_genotype_counts[1] += count_diff;

    for (size_t i=0; i<=2; ++i)
        if (target_genotype_counts[i] > genotype_bins[i].size())
            throw runtime_error("Error: not enough genotypes.");

    cout << "target counts: " << target_genotype_counts[0] << " " << 
        target_genotype_counts[1] <<  " " << target_genotype_counts[2] << endl;

    const Organisms& organisms = population.organisms();
    Organisms organisms_subset;

    vector<size_t> organism_indices;

    for (size_t g=0; g<=2; ++g)
    {
        set<size_t>::const_iterator it=genotype_bins[g].begin();
        for (size_t i=0; i!=target_genotype_counts[g]; ++i, ++it)
        {
            organisms_subset.push_back(organisms[*it]);
            organism_indices.push_back(*it);
        }
    }

    Population_Organisms population_subset(organisms_subset);

    cout << "Writing focal subset population to " << config.filename_subset_population << endl;
    ofstream os(config.filename_subset_population.c_str());
    os << population_subset;
    os.close();

    cout << "Reading haplotypes from " << config.filename_full_ms << endl;

    ifstream is_ms(config.filename_full_ms.c_str());
    vector<string> header_lines;
    vector<string> haplotypes;
    while (is_ms)
    {
        string buffer;
        getline(is_ms, buffer);
        if (!is_ms) break;
        if (buffer.empty() || buffer[0]=='#' || buffer.substr(0,9) == "segsites:")
        {
            header_lines.push_back(buffer);
            continue;
        }
        haplotypes.push_back(buffer);
    }

    cout << haplotypes.size() << " haplotypes read.\n";

    cout << "Writing focal subset haplotypes to " << config.filename_subset_ms << endl;
    ofstream os_ms(config.filename_subset_ms.c_str());
    copy(header_lines.begin(), header_lines.end(), ostream_iterator<string>(os_ms, "\n"));
    for (vector<size_t>::const_iterator it=organism_indices.begin(); it!=organism_indices.end(); ++it)
    {
        os_ms << haplotypes.at(*it * 2) << endl;
        os_ms << haplotypes.at(*it * 2 + 1) << endl;
    }
    os_ms.close();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc < 2)
        {
            cerr << "Usage: forqs_focal_subset <config>\n\n"
                 << "Darren Kessner\n"
                 << "John Novembre Lab, UCLA\n\n";

            Config config;
            config.filename_full_population = "full.pop";
            config.filename_full_ms = "full.ms";
            config.filename_original_snps = "original_snps.ms";
            config.filename_subset_population = "subset.pop";
            config.filename_subset_ms = "subset.ms";
            config.population_size = 100;
            config.focal_locus = LocusPtr(new Locus("dummy", 0, 123456));
            config.focal_allele_frequency = .3;
            cout << "#\n"
                 << "# example forqs_focal_subset config file\n" 
                 << "#\n\n"
                 << config << endl;
            return 1;
        }

        string config_filename = argv[1];

        Config config(config_filename);
        cout << "forqs_focal_subset\n" << config << endl;
        config.validate();

        create_subset_population(config);

        return 0;
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


